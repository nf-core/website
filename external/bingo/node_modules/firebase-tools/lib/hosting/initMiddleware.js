"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const request = require("request");
const logger = require("../logger");
const utils = require("../utils");
const SDK_PATH_REGEXP = /^\/__\/firebase\/([^/]+)\/([^/]+)$/;
function initMiddleware(init) {
    return (req, res, next) => {
        const match = RegExp(SDK_PATH_REGEXP).exec(req.url);
        if (match) {
            const version = match[1];
            const sdkName = match[2];
            const url = `https://www.gstatic.com/firebasejs/${version}/${sdkName}`;
            const preq = request(url)
                .on("response", (pres) => {
                if (pres.statusCode === 404) {
                    return next();
                }
                return preq.pipe(res);
            })
                .on("error", (e) => {
                utils.logLabeledWarning("hosting", `Could not load Firebase SDK ${sdkName} v${version}, check your internet connection.`);
                logger.debug(e);
            });
        }
        else if (req.url === "/__/firebase/init.js") {
            res.setHeader("Content-Type", "application/javascript");
            res.end(init.js);
        }
        else if (req.url === "/__/firebase/init.json") {
            res.setHeader("Content-Type", "application/json");
            res.end(init.json);
        }
        else {
            next();
        }
    };
}
exports.initMiddleware = initMiddleware;
