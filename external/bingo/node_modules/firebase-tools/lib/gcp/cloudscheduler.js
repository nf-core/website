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
const _ = require("lodash");
const api = require("../api");
const error_1 = require("../error");
const utils_1 = require("../utils");
const VERSION = "v1beta1";
const DEFAULT_TIME_ZONE = "America/Los_Angeles";
function createJob(job) {
    const strippedName = job.name.substring(0, job.name.lastIndexOf("/"));
    return api.request("POST", `/${VERSION}/${strippedName}`, {
        auth: true,
        origin: api.cloudschedulerOrigin,
        data: Object.assign({ timeZone: DEFAULT_TIME_ZONE }, job),
    });
}
exports.createJob = createJob;
function deleteJob(name) {
    return api.request("DELETE", `/${VERSION}/${name}`, {
        auth: true,
        origin: api.cloudschedulerOrigin,
    });
}
exports.deleteJob = deleteJob;
function getJob(name) {
    return api.request("GET", `/${VERSION}/${name}`, {
        auth: true,
        origin: api.cloudschedulerOrigin,
        resolveOnHTTPError: true,
    });
}
exports.getJob = getJob;
function updateJob(job) {
    return api.request("PATCH", `/${VERSION}/${job.name}`, {
        auth: true,
        origin: api.cloudschedulerOrigin,
        data: Object.assign({ timeZone: DEFAULT_TIME_ZONE }, job),
    });
}
exports.updateJob = updateJob;
function createOrReplaceJob(job) {
    return __awaiter(this, void 0, void 0, function* () {
        const jobName = job.name.split("/").pop();
        const existingJob = yield getJob(job.name);
        if (existingJob.status === 404) {
            let newJob;
            try {
                newJob = yield createJob(job);
            }
            catch (err) {
                if (_.get(err, "context.response.statusCode") === 404) {
                    return Promise.reject(new error_1.FirebaseError(`Cloud resource location is not set for this project but scheduled functions requires it. ` +
                        `Please see this documentation for more details: https://firebase.google.com/docs/projects/locations.`));
                }
            }
            utils_1.logLabeledSuccess("functions", `created scheduler job ${jobName}`);
            return newJob;
        }
        if (!job.timeZone) {
            job.timeZone = DEFAULT_TIME_ZONE;
        }
        if (isIdentical(existingJob.body, job)) {
            utils_1.logLabeledBullet("functions", `scheduler job ${jobName} is up to date, no changes required`);
            return;
        }
        const updatedJob = yield updateJob(job);
        utils_1.logLabeledBullet("functions", `updated scheduler job ${jobName}`);
        return updatedJob;
    });
}
exports.createOrReplaceJob = createOrReplaceJob;
function isIdentical(job, otherJob) {
    return (job &&
        otherJob &&
        job.schedule === otherJob.schedule &&
        job.timeZone === otherJob.timeZone &&
        _.isEqual(job.retryConfig, otherJob.retryConfig));
}
