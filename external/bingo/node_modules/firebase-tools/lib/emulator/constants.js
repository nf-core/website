"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const url = require("url");
const types_1 = require("./types");
const DEFAULT_PORTS = {
    gui: 4000,
    hub: 4400,
    hosting: 5000,
    functions: 5001,
    firestore: 8080,
    database: 9000,
    pubsub: 8085,
};
const DEFAULT_HOST = "localhost";
class Constants {
    static getServiceName(service) {
        switch (service) {
            case this.SERVICE_FIRESTORE:
                return "firestore";
            case this.SERVICE_REALTIME_DATABASE:
                return "database";
            case this.SERVICE_PUBSUB:
                return "pubsub";
            case this.SERVICE_ANALYTICS:
                return "analytics";
            case this.SERVICE_AUTH:
                return "auth";
            case this.SERVICE_CRASHLYTICS:
                return "crashlytics";
            case this.SERVICE_REMOTE_CONFIG:
                return "remote config";
            case this.SERVICE_STORAGE:
                return "storage";
            case this.SERVICE_TEST_LAB:
                return "test lab";
            default:
                return service;
        }
    }
    static getDefaultHost(emulator) {
        return DEFAULT_HOST;
    }
    static getDefaultPort(emulator) {
        return DEFAULT_PORTS[emulator];
    }
    static getHostKey(emulator) {
        return `emulators.${emulator.toString()}.host`;
    }
    static getPortKey(emulator) {
        return `emulators.${emulator.toString()}.port`;
    }
    static getAddress(emulator, options) {
        const hostVal = options.config.get(this.getHostKey(emulator), DEFAULT_HOST);
        const portVal = options.config.get(this.getPortKey(emulator), this.getDefaultPort(emulator));
        const host = this.normalizeHost(hostVal);
        const port = parseInt(portVal, 10);
        return { host, port };
    }
    static description(name) {
        if (name === types_1.Emulators.HUB) {
            return "emulator hub";
        }
        else if (name === types_1.Emulators.GUI) {
            return "emulator GUI";
        }
        else {
            return `${name} emulator`;
        }
    }
    static normalizeHost(host) {
        let normalized = host;
        if (!normalized.startsWith("http")) {
            normalized = `http://${normalized}`;
        }
        const u = url.parse(normalized);
        return u.hostname || DEFAULT_HOST;
    }
}
exports.Constants = Constants;
Constants.DEFAULT_DATABASE_EMULATOR_NAMESPACE = "fake-server";
Constants.SERVICE_FIRESTORE = "firestore.googleapis.com";
Constants.SERVICE_REALTIME_DATABASE = "firebaseio.com";
Constants.SERVICE_PUBSUB = "pubsub.googleapis.com";
Constants.SERVICE_ANALYTICS = "app-measurement.com";
Constants.SERVICE_AUTH = "firebaseauth.googleapis.com";
Constants.SERVICE_CRASHLYTICS = "fabric.io";
Constants.SERVICE_REMOTE_CONFIG = "firebaseremoteconfig.googleapis.com";
Constants.SERVICE_STORAGE = "storage.googleapis.com";
Constants.SERVICE_TEST_LAB = "testing.googleapis.com";
