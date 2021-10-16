"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const constants_1 = require("../../emulator/constants");
const emulatorLogger_1 = require("../../emulator/emulatorLogger");
function functionResourceToEmulatedTriggerDefintion(resource) {
    const etd = {
        name: resource.name,
        entryPoint: resource.name,
    };
    const properties = _.get(resource, "properties", {});
    if (properties.timeout) {
        etd.timeout = properties.timeout;
    }
    if (properties.location) {
        etd.regions = [properties.location];
    }
    if (properties.availableMemoryMb) {
        etd.availableMemoryMb = properties.availableMemoryMb;
    }
    if (properties.httpsTrigger) {
        etd.httpsTrigger = properties.httpsTrigger;
    }
    else if (properties.eventTrigger) {
        properties.eventTrigger.service = getServiceFromEventType(properties.eventTrigger.eventType);
        etd.eventTrigger = properties.eventTrigger;
    }
    else {
        emulatorLogger_1.EmulatorLogger.log("WARN", `Function '${resource.name} is missing a trigger in extension.yaml. Please add one, as triggers defined in code are ignored.`);
    }
    return etd;
}
exports.functionResourceToEmulatedTriggerDefintion = functionResourceToEmulatedTriggerDefintion;
function getServiceFromEventType(eventType) {
    if (eventType.includes("firestore")) {
        return constants_1.Constants.SERVICE_FIRESTORE;
    }
    if (eventType.includes("database")) {
        return constants_1.Constants.SERVICE_REALTIME_DATABASE;
    }
    if (eventType.includes("pubsub")) {
        return constants_1.Constants.SERVICE_PUBSUB;
    }
    if (eventType.includes("analytics")) {
        return constants_1.Constants.SERVICE_ANALYTICS;
    }
    if (eventType.includes("auth")) {
        return constants_1.Constants.SERVICE_AUTH;
    }
    if (eventType.includes("crashlytics")) {
        return constants_1.Constants.SERVICE_CRASHLYTICS;
    }
    if (eventType.includes("remoteconfig")) {
        return constants_1.Constants.SERVICE_REMOTE_CONFIG;
    }
    if (eventType.includes("storage")) {
        return constants_1.Constants.SERVICE_STORAGE;
    }
    if (eventType.includes("testing")) {
        return constants_1.Constants.SERVICE_TEST_LAB;
    }
    return "";
}
