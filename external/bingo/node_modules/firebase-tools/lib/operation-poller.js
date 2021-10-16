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
const queue_1 = require("./throttler/queue");
const api = require("./api");
const error_1 = require("./error");
const DEFAULT_INITIAL_BACKOFF_DELAY_MILLIS = 250;
const DEFAULT_MASTER_TIMEOUT_MILLIS = 30000;
class OperationPoller {
    poll(options) {
        return __awaiter(this, void 0, void 0, function* () {
            const queue = new queue_1.Queue({
                name: options.pollerName || "LRO Poller",
                concurrency: 1,
                retries: Number.MAX_SAFE_INTEGER,
                backoff: options.backoff || DEFAULT_INITIAL_BACKOFF_DELAY_MILLIS,
            });
            const masterTimeout = options.masterTimeout || DEFAULT_MASTER_TIMEOUT_MILLIS;
            const { response, error } = yield queue.run(this.getPollingTask(options), masterTimeout);
            queue.close();
            if (error) {
                throw error instanceof error_1.FirebaseError
                    ? error
                    : new error_1.FirebaseError(error.message, { status: error.code });
            }
            return response;
        });
    }
    getPollingTask(options) {
        const { apiOrigin, apiVersion, operationResourceName } = options;
        const requestOptions = { auth: true, origin: apiOrigin };
        return () => __awaiter(this, void 0, void 0, function* () {
            let apiResponse;
            try {
                apiResponse = yield api.request("GET", `/${apiVersion}/${operationResourceName}`, requestOptions);
            }
            catch (err) {
                if (err.status === 500 || err.status === 503) {
                    throw err;
                }
                return { error: err };
            }
            if (!apiResponse.body || !apiResponse.body.done) {
                throw new Error("Polling incomplete, should trigger retry with backoff");
            }
            return apiResponse.body;
        });
    }
}
exports.OperationPoller = OperationPoller;
function pollOperation(options) {
    return new OperationPoller().poll(options);
}
exports.pollOperation = pollOperation;
