"use strict";
var clc = require("cli-color");
var ProgressBar = require("progress");
var api = require("../api");
var firestore = require("../gcp/firestore");
var { FirebaseError } = require("../error");
var logger = require("../logger");
var utils = require("../utils");
var MIN_ID = "__id-9223372036854775808__";
function FirestoreDelete(project, path, options) {
    this.project = project;
    this.path = path || "";
    this.recursive = Boolean(options.recursive);
    this.shallow = Boolean(options.shallow);
    this.allCollections = Boolean(options.allCollections);
    this.path = this.path.replace(/(^\/+|\/+$)/g, "");
    this.allDescendants = this.recursive;
    this.root = "projects/" + project + "/databases/(default)/documents";
    var segments = this.path.split("/");
    this.isDocumentPath = segments.length % 2 === 0;
    this.isCollectionPath = !this.isDocumentPath;
    this.parent = this.root;
    if (this.isCollectionPath) {
        segments.pop();
    }
    if (segments.length > 0) {
        this.parent += "/" + segments.join("/");
    }
    if (!options.allCollections) {
        this._validateOptions();
    }
}
FirestoreDelete.prototype._validateOptions = function () {
    if (this.recursive && this.shallow) {
        throw new FirebaseError("Cannot pass recursive and shallow options together.");
    }
    if (this.isCollectionPath && !this.recursive && !this.shallow) {
        throw new FirebaseError("Must pass recursive or shallow option when deleting a collection.");
    }
    var pieces = this.path.split("/");
    if (pieces.length === 0) {
        throw new FirebaseError("Path length must be greater than zero.");
    }
    var hasEmptySegment = pieces.some(function (piece) {
        return piece.length === 0;
    });
    if (hasEmptySegment) {
        throw new FirebaseError("Path must not have any empty segments.");
    }
};
FirestoreDelete.prototype._collectionDescendantsQuery = function (allDescendants, batchSize, startAfter) {
    var nullChar = String.fromCharCode(0);
    var startAt = this.root + "/" + this.path + "/" + MIN_ID;
    var endAt = this.root + "/" + this.path + nullChar + "/" + MIN_ID;
    var where = {
        compositeFilter: {
            op: "AND",
            filters: [
                {
                    fieldFilter: {
                        field: {
                            fieldPath: "__name__",
                        },
                        op: "GREATER_THAN_OR_EQUAL",
                        value: {
                            referenceValue: startAt,
                        },
                    },
                },
                {
                    fieldFilter: {
                        field: {
                            fieldPath: "__name__",
                        },
                        op: "LESS_THAN",
                        value: {
                            referenceValue: endAt,
                        },
                    },
                },
            ],
        },
    };
    var query = {
        structuredQuery: {
            where: where,
            limit: batchSize,
            from: [
                {
                    allDescendants: allDescendants,
                },
            ],
            select: {
                fields: [{ fieldPath: "__name__" }],
            },
            orderBy: [{ field: { fieldPath: "__name__" } }],
        },
    };
    if (startAfter) {
        query.structuredQuery.startAt = {
            values: [{ referenceValue: startAfter }],
            before: false,
        };
    }
    return query;
};
FirestoreDelete.prototype._docDescendantsQuery = function (allDescendants, batchSize, startAfter) {
    var query = {
        structuredQuery: {
            limit: batchSize,
            from: [
                {
                    allDescendants: allDescendants,
                },
            ],
            select: {
                fields: [{ fieldPath: "__name__" }],
            },
            orderBy: [{ field: { fieldPath: "__name__" } }],
        },
    };
    if (startAfter) {
        query.structuredQuery.startAt = {
            values: [{ referenceValue: startAfter }],
            before: false,
        };
    }
    return query;
};
FirestoreDelete.prototype._getDescendantBatch = function (allDescendants, batchSize, startAfter) {
    var url = this.parent + ":runQuery";
    var body;
    if (this.isDocumentPath) {
        body = this._docDescendantsQuery(allDescendants, batchSize, startAfter);
    }
    else {
        body = this._collectionDescendantsQuery(allDescendants, batchSize, startAfter);
    }
    return api
        .request("POST", "/v1beta1/" + url, {
        auth: true,
        data: body,
        origin: api.firestoreOrigin,
    })
        .then(function (res) {
        return res.body
            .filter(function (x) {
            return x.document;
        })
            .map(function (x) {
            return x.document;
        });
    });
};
FirestoreDelete.progressBar = new ProgressBar("Deleted :current docs (:rate docs/s)\n", {
    total: Number.MAX_SAFE_INTEGER,
});
FirestoreDelete.prototype._recursiveBatchDelete = function () {
    var self = this;
    var readBatchSize = 7500;
    var deleteBatchSize = 250;
    var maxPendingDeletes = 15;
    var maxQueueSize = deleteBatchSize * maxPendingDeletes * 2;
    var queue = [];
    var numPendingDeletes = 0;
    var pagesRemaining = true;
    var pageIncoming = false;
    var lastDocName;
    var failures = [];
    var retried = {};
    var queueLoop = function () {
        if (queue.length == 0 && numPendingDeletes == 0 && !pagesRemaining) {
            return true;
        }
        if (failures.length > 0) {
            logger.debug("Found " + failures.length + " failed deletes, failing.");
            return true;
        }
        if (queue.length <= maxQueueSize && pagesRemaining && !pageIncoming) {
            pageIncoming = true;
            self
                ._getDescendantBatch(self.allDescendants, readBatchSize, lastDocName)
                .then(function (docs) {
                pageIncoming = false;
                if (docs.length == 0) {
                    pagesRemaining = false;
                    return;
                }
                queue = queue.concat(docs);
                lastDocName = docs[docs.length - 1].name;
            })
                .catch(function (e) {
                logger.debug("Failed to fetch page after " + lastDocName, e);
                pageIncoming = false;
            });
        }
        if (numPendingDeletes > maxPendingDeletes) {
            return false;
        }
        if (queue.length == 0) {
            return false;
        }
        var toDelete = [];
        var numToDelete = Math.min(deleteBatchSize, queue.length);
        for (var i = 0; i < numToDelete; i++) {
            toDelete.push(queue.shift());
        }
        numPendingDeletes++;
        firestore
            .deleteDocuments(self.project, toDelete)
            .then(function (numDeleted) {
            FirestoreDelete.progressBar.tick(numDeleted);
            numPendingDeletes--;
        })
            .catch(function (e) {
            if (e.status >= 500 && e.status < 600) {
                logger.debug("Server error deleting doc batch", e);
                toDelete.forEach(function (doc) {
                    if (retried[doc.name]) {
                        logger.debug("Failed to delete doc " + doc.name + " multiple times.");
                        failures.push(doc.name);
                    }
                    else {
                        retried[doc.name] = true;
                        queue.push(doc);
                    }
                });
            }
            else {
                logger.debug("Fatal error deleting docs ", e);
                failures = failures.concat(toDelete);
            }
            numPendingDeletes--;
        });
        return false;
    };
    return new Promise(function (resolve, reject) {
        var intervalId = setInterval(function () {
            if (queueLoop()) {
                clearInterval(intervalId);
                if (failures.length == 0) {
                    resolve();
                }
                else {
                    reject("Failed to delete documents " + failures);
                }
            }
        }, 0);
    });
};
FirestoreDelete.prototype._deletePath = function () {
    var self = this;
    var initialDelete;
    if (this.isDocumentPath) {
        var doc = { name: this.root + "/" + this.path };
        initialDelete = firestore.deleteDocument(doc).catch(function (err) {
            logger.debug("deletePath:initialDelete:error", err);
            if (self.allDescendants) {
                return Promise.resolve();
            }
            return utils.reject("Unable to delete " + clc.cyan(this.path));
        });
    }
    else {
        initialDelete = Promise.resolve();
    }
    return initialDelete.then(function () {
        return self._recursiveBatchDelete();
    });
};
FirestoreDelete.prototype.deleteDatabase = function () {
    var self = this;
    return firestore
        .listCollectionIds(this.project)
        .catch(function (err) {
        logger.debug("deleteDatabase:listCollectionIds:error", err);
        return utils.reject("Unable to list collection IDs");
    })
        .then(function (collectionIds) {
        var promises = [];
        logger.info("Deleting the following collections: " + clc.cyan(collectionIds.join(", ")));
        for (var i = 0; i < collectionIds.length; i++) {
            var collectionId = collectionIds[i];
            var deleteOp = new FirestoreDelete(self.project, collectionId, {
                recursive: true,
            });
            promises.push(deleteOp.execute());
        }
        return Promise.all(promises);
    });
};
FirestoreDelete.prototype.checkHasChildren = function () {
    return this._getDescendantBatch(true, 1).then(function (docs) {
        return docs.length > 0;
    });
};
FirestoreDelete.prototype.execute = function () {
    var verifyRecurseSafe;
    if (this.isDocumentPath && !this.recursive && !this.shallow) {
        verifyRecurseSafe = this.checkHasChildren().then(function (multiple) {
            if (multiple) {
                return utils.reject("Document has children, must specify -r or --shallow.", { exit: 1 });
            }
        });
    }
    else {
        verifyRecurseSafe = Promise.resolve();
    }
    var self = this;
    return verifyRecurseSafe.then(function () {
        return self._deletePath();
    });
};
module.exports = FirestoreDelete;
