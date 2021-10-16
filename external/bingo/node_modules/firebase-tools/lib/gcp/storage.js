"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var path = require("path");
var api = require("../api");
var logger = require("../logger");
var { FirebaseError } = require("../error");
function _getDefaultBucket(projectId) {
    return api
        .request("GET", "/v1/apps/" + projectId, {
        auth: true,
        origin: api.appengineOrigin,
    })
        .then(function (resp) {
        if (resp.body.defaultBucket === "undefined") {
            logger.debug("Default storage bucket is undefined.");
            return Promise.reject(new FirebaseError("Your project is being set up. Please wait a minute before deploying again."));
        }
        return Promise.resolve(resp.body.defaultBucket);
    }, function (err) {
        logger.info("\n\nThere was an issue deploying your functions. Verify that your project has a Google App Engine instance setup at https://console.cloud.google.com/appengine and try again. If this issue persists, please contact support.");
        return Promise.reject(err);
    });
}
function _uploadSource(source, uploadUrl) {
    return api.request("PUT", uploadUrl, {
        data: source.stream,
        headers: {
            "Content-Type": "application/zip",
            "x-goog-content-length-range": "0,104857600",
        },
        json: false,
        origin: api.storageOrigin,
        logOptions: { skipRequestBody: true },
    });
}
function _uploadObject(source, bucketName) {
    return __awaiter(this, void 0, void 0, function* () {
        if (path.extname(source.file) !== ".zip") {
            throw new FirebaseError(`Expected a file name ending in .zip, got ${source.file}`);
        }
        const location = `/v0/b/${bucketName}/o/${path.basename(source.file)}`;
        yield api.request("POST", location, {
            auth: true,
            data: source.stream,
            headers: {
                "Content-Type": "application/zip",
                "x-goog-content-length-range": "0,123289600",
            },
            json: false,
            origin: api.firebaseStorageOrigin,
            logOptions: { skipRequestBody: true },
        });
        return location;
    });
}
function _deleteObject(location) {
    return api.request("DELETE", location, {
        auth: true,
        origin: api.firebaseStorageOrigin,
    });
}
module.exports = {
    getDefaultBucket: _getDefaultBucket,
    deleteObject: _deleteObject,
    upload: _uploadSource,
    uploadObject: _uploadObject,
};
