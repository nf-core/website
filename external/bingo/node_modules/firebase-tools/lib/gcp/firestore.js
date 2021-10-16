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
Object.defineProperty(exports, "__esModule", { value: true });
const api = require("../api");
const _API_ROOT = "/v1beta1/";
function listCollectionIds(project) {
    return __awaiter(this, void 0, void 0, function* () {
        const url = _API_ROOT + "projects/" + project + "/databases/(default)/documents:listCollectionIds";
        return api
            .request("POST", url, {
            auth: true,
            origin: api.firestoreOrigin,
            data: {
                pageSize: 2147483647,
            },
        })
            .then((res) => {
            return res.body.collectionIds || [];
        });
    });
}
exports.listCollectionIds = listCollectionIds;
function deleteDocument(doc) {
    return __awaiter(this, void 0, void 0, function* () {
        return api.request("DELETE", _API_ROOT + doc.name, {
            auth: true,
            origin: api.firestoreOrigin,
        });
    });
}
exports.deleteDocument = deleteDocument;
function deleteDocuments(project, docs) {
    return __awaiter(this, void 0, void 0, function* () {
        const url = _API_ROOT + "projects/" + project + "/databases/(default)/documents:commit";
        const writes = docs.map((doc) => {
            return { delete: doc.name };
        });
        const body = { writes };
        try {
            const res = yield api.request("POST", url, {
                auth: true,
                data: body,
                origin: api.firestoreOrigin,
            });
            return res.body.writeResults.length;
        }
        catch (err) {
            if (err.status === 400 &&
                err.message.indexOf("Transaction too big") !== -1 &&
                docs.length > 2) {
                const a = yield deleteDocuments(project, docs.slice(0, docs.length / 2));
                const b = yield deleteDocuments(project, docs.slice(docs.length / 2));
                return a + b;
            }
            throw err;
        }
    });
}
exports.deleteDocuments = deleteDocuments;
