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
const uuid = require("uuid");
const pubsub_1 = require("@google-cloud/pubsub");
const api = require("../api");
const downloadableEmulators = require("./downloadableEmulators");
const emulatorLogger_1 = require("./emulatorLogger");
const types_1 = require("../emulator/types");
const constants_1 = require("./constants");
const error_1 = require("../error");
const registry_1 = require("./registry");
class PubsubEmulator {
    constructor(args) {
        this.args = args;
        const { host, port } = this.getInfo();
        this.pubsub = new pubsub_1.PubSub({
            apiEndpoint: `${host}:${port}`,
            projectId: this.args.projectId,
        });
        this.triggers = new Map();
        this.subscriptions = new Map();
    }
    start() {
        return __awaiter(this, void 0, void 0, function* () {
            return downloadableEmulators.start(types_1.Emulators.PUBSUB, this.args);
        });
    }
    connect() {
        return __awaiter(this, void 0, void 0, function* () {
            return;
        });
    }
    stop() {
        return __awaiter(this, void 0, void 0, function* () {
            yield downloadableEmulators.stop(types_1.Emulators.PUBSUB);
        });
    }
    getInfo() {
        const host = this.args.host || constants_1.Constants.getDefaultHost(types_1.Emulators.PUBSUB);
        const port = this.args.port || constants_1.Constants.getDefaultPort(types_1.Emulators.PUBSUB);
        return {
            host,
            port,
        };
    }
    getName() {
        return types_1.Emulators.PUBSUB;
    }
    addTrigger(topicName, trigger) {
        return __awaiter(this, void 0, void 0, function* () {
            emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `addTrigger(${topicName}, ${trigger})`);
            const topicTriggers = this.triggers.get(topicName) || new Set();
            if (topicTriggers.has(topicName) && this.subscriptions.has(topicName)) {
                emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", "Trigger already exists");
                return;
            }
            const topic = this.pubsub.topic(topicName);
            try {
                emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `Creating topic: ${topicName}`);
                yield topic.create();
            }
            catch (e) {
                if (e && e.code === 6) {
                    emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `Topic ${topicName} exists`);
                }
                else {
                    throw new error_1.FirebaseError(`Could not create topic ${topicName}`, { original: e });
                }
            }
            const subName = `emulator-sub-${topicName}`;
            let sub;
            try {
                emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `Creating sub for topic: ${topicName}`);
                [sub] = yield topic.createSubscription(subName);
            }
            catch (e) {
                if (e && e.code === 6) {
                    emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `Sub for ${topicName} exists`);
                    sub = topic.subscription(`emulator-sub-${topicName}`);
                }
                else {
                    throw new error_1.FirebaseError(`Could not create sub ${subName}`, { original: e });
                }
            }
            sub.on("message", (message) => {
                this.onMessage(topicName, message);
            });
            topicTriggers.add(trigger);
            this.triggers.set(topicName, topicTriggers);
            this.subscriptions.set(topicName, sub);
        });
    }
    onMessage(topicName, message) {
        return __awaiter(this, void 0, void 0, function* () {
            emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `onMessage(${topicName}, ${message.id})`);
            const topicTriggers = this.triggers.get(topicName);
            if (!topicTriggers || topicTriggers.size === 0) {
                throw new error_1.FirebaseError(`No trigger for topic: ${topicName}`);
            }
            const functionsEmu = registry_1.EmulatorRegistry.get(types_1.Emulators.FUNCTIONS);
            if (!functionsEmu) {
                throw new error_1.FirebaseError(`Attempted to execute pubsub trigger for topic ${topicName} but could not find Functions emulator`);
            }
            const functionsPort = functionsEmu.getInfo().port;
            const functionsHost = functionsEmu.getInfo().host;
            emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `Executing ${topicTriggers.size} matching triggers (${JSON.stringify(Array.from(topicTriggers))})`);
            let remaining = topicTriggers.size;
            for (const trigger of topicTriggers) {
                const body = {
                    context: {
                        eventId: uuid.v4(),
                        resource: {
                            service: "pubsub.googleapis.com",
                            name: `projects/${this.args.projectId}/topics/${topicName}`,
                        },
                        eventType: "google.pubsub.topic.publish",
                        timestamp: message.publishTime.toISOString(),
                    },
                    data: {
                        data: message.data,
                        attributes: message.attributes,
                    },
                };
                try {
                    yield api.request("POST", `/functions/projects/${this.args.projectId}/triggers/${trigger}`, {
                        origin: `http://${functionsHost}:${functionsPort}`,
                        data: body,
                    });
                }
                catch (e) {
                    emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", e);
                }
                remaining--;
                if (remaining <= 0) {
                    emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "pubsub", `Acking message ${message.id}`);
                    message.ack();
                }
            }
        });
    }
}
exports.PubsubEmulator = PubsubEmulator;
