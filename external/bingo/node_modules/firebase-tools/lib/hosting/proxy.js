"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const lodash_1 = require("lodash");
const request = require("request");
const logger = require("../logger");
const REQUIRED_VARY_VALUES = ["Accept-Encoding", "Authorization", "Cookie"];
function makeVary(vary) {
    if (!vary) {
        return "Accept-Encoding, Authorization, Cookie";
    }
    const varies = vary.split(/, ?/).map((v) => {
        return v
            .split("-")
            .map((part) => lodash_1.capitalize(part))
            .join("-");
    });
    REQUIRED_VARY_VALUES.forEach((requiredVary) => {
        if (!lodash_1.includes(varies, requiredVary)) {
            varies.push(requiredVary);
        }
    });
    return varies.join(", ");
}
function proxyRequestHandler(url, rewriteIdentifier) {
    return (req, res, next) => {
        logger.info(`[hosting] Rewriting ${req.url} to ${url} for ${rewriteIdentifier}`);
        const cookie = req.headers.cookie || "";
        const sessionCookie = cookie.split(/; ?/).find((c) => {
            return c.trim().indexOf("__session=") === 0;
        });
        const proxied = request({
            method: req.method,
            qs: req.query,
            url: url + req.url,
            headers: {
                "X-Forwarded-Host": req.headers.host,
                "X-Original-Url": req.url,
                Pragma: "no-cache",
                "Cache-Control": "no-cache, no-store",
                Cookie: sessionCookie,
            },
            followRedirect: false,
            timeout: 60000,
        });
        req.pipe(proxied);
        proxied.on("error", (err) => {
            if (err.code === "ETIMEDOUT" || err.code === "ESOCKETTIMEDOUT") {
                res.statusCode = 504;
                return res.end("Timed out waiting for function to respond.");
            }
            res.statusCode = 500;
            res.end(`An internal error occurred while proxying for ${rewriteIdentifier}`);
        });
        proxied.on("response", (response) => {
            if (response.statusCode === 404) {
                const cascade = response.headers["x-cascade"];
                if (cascade && cascade.toUpperCase() === "PASS") {
                    return next();
                }
            }
            if (!response.headers["cache-control"]) {
                response.headers["cache-control"] = "private";
            }
            if (response.headers["cache-control"].indexOf("private") < 0) {
                delete response.headers["set-cookie"];
            }
            response.headers.vary = makeVary(response.headers.vary);
            proxied.pipe(res);
        });
    };
}
exports.proxyRequestHandler = proxyRequestHandler;
function errorRequestHandler(error) {
    return (req, res, next) => {
        res.statusCode = 500;
        const out = `A problem occurred while trying to handle a proxied rewrite: ${error}`;
        logger.error(out);
        res.end(out);
    };
}
exports.errorRequestHandler = errorRequestHandler;
