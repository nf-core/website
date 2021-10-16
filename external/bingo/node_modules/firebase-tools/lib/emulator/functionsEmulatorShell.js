"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const uuid = require("uuid");
const functionsEmulator_1 = require("./functionsEmulator");
const functionsEmulatorShared_1 = require("./functionsEmulatorShared");
const utils = require("../utils");
const logger = require("../logger");
const error_1 = require("../error");
class FunctionsEmulatorShell {
    constructor(emu) {
        this.emu = emu;
        this.urls = {};
        this.triggers = emu.getTriggers();
        this.emulatedFunctions = this.triggers.map((t) => t.name);
        const entryPoints = this.triggers.map((t) => t.entryPoint);
        utils.logLabeledBullet("functions", `Loaded functions: ${entryPoints.join(", ")}`);
        for (const trigger of this.triggers) {
            const name = trigger.name;
            if (trigger.httpsTrigger) {
                this.urls[name] = functionsEmulator_1.FunctionsEmulator.getHttpFunctionUrl(this.emu.getInfo().host, this.emu.getInfo().port, this.emu.getProjectId(), name, functionsEmulatorShared_1.getFunctionRegion(trigger));
            }
        }
    }
    call(name, data, opts) {
        const trigger = this.getTrigger(name);
        logger.debug(`shell:${name}: trigger=${JSON.stringify(trigger)}`);
        logger.debug(`shell:${name}: opts=${JSON.stringify(opts)}, data=${JSON.stringify(data)}`);
        if (!trigger.eventTrigger) {
            throw new error_1.FirebaseError(`Function ${name} is not a background function`);
        }
        const eventType = trigger.eventTrigger.eventType;
        let resource = opts.resource;
        if (typeof resource === "object" && resource.name) {
            resource = resource.name;
        }
        const proto = {
            eventId: uuid.v4(),
            timestamp: new Date().toISOString(),
            eventType,
            resource,
            params: opts.params,
            auth: opts.auth,
            data,
        };
        this.emu.startFunctionRuntime(name, functionsEmulatorShared_1.EmulatedTriggerType.BACKGROUND, proto);
    }
    getTrigger(name) {
        const result = this.triggers.find((trigger) => {
            return trigger.name === name;
        });
        if (!result) {
            throw new error_1.FirebaseError(`Could not find trigger ${name}`);
        }
        return result;
    }
}
exports.FunctionsEmulatorShell = FunctionsEmulatorShell;
