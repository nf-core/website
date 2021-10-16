"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const api = require("../api");
const VERSION = "v1";
function createTopic(name) {
    return api.request("PUT", `/${VERSION}/${name}`, {
        auth: true,
        origin: api.pubsubOrigin,
        data: { labels: { deployment: "firebase-schedule" } },
    });
}
exports.createTopic = createTopic;
function deleteTopic(name) {
    return api.request("DELETE", `/${VERSION}/${name}`, {
        auth: true,
        origin: api.pubsubOrigin,
    });
}
exports.deleteTopic = deleteTopic;
