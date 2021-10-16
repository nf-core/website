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
const ensureApiEnabled_1 = require("../../ensureApiEnabled");
const functionsDeployHelper_1 = require("../../functionsDeployHelper");
const cloudscheduler_1 = require("../../gcp/cloudscheduler");
const pubsub_1 = require("../../gcp/pubsub");
function createOrUpdateSchedulesAndTopics(projectId, triggers, existingScheduledFunctions, appEngineLocation) {
    return __awaiter(this, void 0, void 0, function* () {
        const triggersWithSchedules = triggers.filter((t) => !!t.schedule);
        let schedulerEnabled = false;
        if (triggersWithSchedules.length) {
            yield Promise.all([
                ensureApiEnabled_1.ensure(projectId, "cloudscheduler.googleapis.com", "scheduler", false),
                ensureApiEnabled_1.ensure(projectId, "pubsub.googleapis.com", "pubsub", false),
            ]);
            schedulerEnabled = true;
        }
        else {
            schedulerEnabled = yield ensureApiEnabled_1.check(projectId, "cloudscheduler.googleapis.com", "scheduler", false);
        }
        for (const trigger of triggers) {
            const [, , , region, , functionName] = trigger.name.split("/");
            const scheduleName = functionsDeployHelper_1.getScheduleName(trigger.name, appEngineLocation);
            const topicName = functionsDeployHelper_1.getTopicName(trigger.name);
            if (!trigger.schedule) {
                if (schedulerEnabled && _.includes(existingScheduledFunctions, trigger.name)) {
                    try {
                        yield cloudscheduler_1.deleteJob(scheduleName);
                    }
                    catch (err) {
                        if (err.context.response.statusCode !== 404) {
                            throw err;
                        }
                    }
                    try {
                        yield pubsub_1.deleteTopic(topicName);
                    }
                    catch (err) {
                        if (err.context.response.statusCode !== 404) {
                            throw err;
                        }
                    }
                }
            }
            else {
                yield cloudscheduler_1.createOrReplaceJob(Object.assign(trigger.schedule, {
                    name: `projects/${projectId}/locations/${appEngineLocation}/jobs/firebase-schedule-${functionName}-${region}`,
                    pubsubTarget: {
                        topicName: `projects/${projectId}/topics/firebase-schedule-${functionName}-${region}`,
                        attributes: {
                            scheduled: "true",
                        },
                    },
                }));
            }
        }
        return;
    });
}
exports.createOrUpdateSchedulesAndTopics = createOrUpdateSchedulesAndTopics;
