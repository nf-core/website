"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const request = require("request");
const utils = require("../utils");
const error_1 = require("../error");
const logger = require("../logger");
const api = require("../api");
class RTDBRemoveRemote {
    constructor(instance) {
        this.instance = instance;
    }
    deletePath(path) {
        return this.patch(path, null, "all data");
    }
    deleteSubPath(path, subPaths) {
        const body = {};
        for (const c of subPaths) {
            body[c] = null;
        }
        return this.patch(path, body, `${subPaths.length} subpaths`);
    }
    patch(path, body, note) {
        const t0 = Date.now();
        return new Promise((resolve, reject) => {
            const url = utils.addSubdomain(api.realtimeOrigin, this.instance) +
                path +
                ".json?print=silent&writeSizeLimit=tiny";
            return api
                .addRequestHeaders({
                url,
                body,
                json: true,
            })
                .then((reqOptionsWithToken) => {
                request.patch(reqOptionsWithToken, (err, res, resBody) => {
                    if (err) {
                        return reject(new error_1.FirebaseError(`Unexpected error while removing data at ${path}`, {
                            exit: 2,
                            original: err,
                        }));
                    }
                    const dt = Date.now() - t0;
                    if (res.statusCode >= 400) {
                        logger.debug(`[database] Failed to remove ${note} at ${path} time: ${dt}ms, will try recursively chunked deletes.`);
                        return resolve(false);
                    }
                    logger.debug(`[database] Sucessfully removed ${note} at ${path} time: ${dt}ms`);
                    return resolve(true);
                });
            });
        });
    }
}
exports.RTDBRemoveRemote = RTDBRemoveRemote;
