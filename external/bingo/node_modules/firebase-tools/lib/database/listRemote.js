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
const request = require("request");
const responseToError = require("../responseToError");
const utils = require("../utils");
const error_1 = require("../error");
const logger = require("../logger");
const api = require("../api");
class RTDBListRemote {
    constructor(instance) {
        this.instance = instance;
    }
    listPath(path, numSubPath, startAfter, timeout) {
        return __awaiter(this, void 0, void 0, function* () {
            const url = `${utils.addSubdomain(api.realtimeOrigin, this.instance)}${path}.json`;
            const params = {
                shallow: true,
                limitToFirst: numSubPath,
            };
            if (startAfter) {
                params.startAfter = startAfter;
            }
            if (timeout) {
                params.timeout = `${timeout}ms`;
            }
            const t0 = Date.now();
            const reqOptionsWithToken = yield api.addRequestHeaders({ url });
            reqOptionsWithToken.qs = params;
            const paths = yield new Promise((resolve, reject) => {
                request.get(reqOptionsWithToken, (err, res, body) => {
                    if (err) {
                        return reject(new error_1.FirebaseError("Unexpected error while listing subtrees", {
                            exit: 2,
                            original: err,
                        }));
                    }
                    else if (res.statusCode >= 400) {
                        return reject(responseToError(res, body));
                    }
                    let data;
                    try {
                        data = JSON.parse(body);
                    }
                    catch (e) {
                        return reject(new error_1.FirebaseError("Malformed JSON response in shallow get ", {
                            exit: 2,
                            original: e,
                        }));
                    }
                    if (data) {
                        return resolve(Object.keys(data));
                    }
                    return resolve([]);
                });
            });
            const dt = Date.now() - t0;
            logger.debug(`[database] sucessfully fetched ${paths.length} path at ${path} ${dt}`);
            return paths;
        });
    }
}
exports.RTDBListRemote = RTDBListRemote;
