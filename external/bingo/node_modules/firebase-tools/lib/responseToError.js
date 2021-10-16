"use strict";
const _ = require("lodash");
const { FirebaseError } = require("./error");
module.exports = function (response, body) {
    if (typeof body === "string" && response.statusCode === 404) {
        body = {
            error: {
                message: "Not Found",
            },
        };
    }
    if (response.statusCode < 400) {
        return null;
    }
    if (typeof body !== "object") {
        try {
            body = JSON.parse(body);
        }
        catch (e) {
            body = {};
        }
    }
    if (!body.error) {
        const errMessage = response.statusCode === 404 ? "Not Found" : "Unknown Error";
        body.error = {
            message: errMessage,
        };
    }
    const message = "HTTP Error: " + response.statusCode + ", " + (body.error.message || body.error);
    let exitCode;
    if (response.statusCode >= 500) {
        exitCode = 2;
    }
    else {
        exitCode = 1;
    }
    _.unset(response, "request.headers");
    return new FirebaseError(message, {
        context: {
            body: body,
            response: response,
        },
        exit: exitCode,
        status: response.statusCode,
    });
};
