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
const types_1 = require("./types");
const functionsEmulatorShared_1 = require("./functionsEmulatorShared");
const events_1 = require("events");
const emulatorLogger_1 = require("./emulatorLogger");
const error_1 = require("../error");
var RuntimeWorkerState;
(function (RuntimeWorkerState) {
    RuntimeWorkerState["IDLE"] = "IDLE";
    RuntimeWorkerState["BUSY"] = "BUSY";
    RuntimeWorkerState["FINISHING"] = "FINISHING";
    RuntimeWorkerState["FINISHED"] = "FINISHED";
})(RuntimeWorkerState = exports.RuntimeWorkerState || (exports.RuntimeWorkerState = {}));
class RuntimeWorker {
    constructor(key, runtime) {
        this.stateEvents = new events_1.EventEmitter();
        this.logListeners = [];
        this._state = RuntimeWorkerState.IDLE;
        this.id = uuid.v4();
        this.key = key;
        this.runtime = runtime;
        this.runtime.events.on("log", (log) => {
            if (log.type === "runtime-status") {
                if (log.data.state === "idle") {
                    if (this.state === RuntimeWorkerState.BUSY) {
                        this.state = RuntimeWorkerState.IDLE;
                    }
                    else if (this.state === RuntimeWorkerState.FINISHING) {
                        this.log(`IDLE --> FINISHING`);
                        this.runtime.shutdown();
                    }
                }
            }
        });
        this.runtime.exit.then(() => {
            this.log("exited");
            this.state = RuntimeWorkerState.FINISHED;
        });
    }
    execute(frb, opts) {
        return __awaiter(this, void 0, void 0, function* () {
            const execFrb = Object.assign({}, frb);
            if (!execFrb.socketPath) {
                execFrb.socketPath = functionsEmulatorShared_1.getTemporarySocketPath(this.runtime.pid, execFrb.cwd);
                this.log(`Assigning socketPath: ${execFrb.socketPath}`);
            }
            const args = { frb: execFrb, opts };
            this.state = RuntimeWorkerState.BUSY;
            this.lastArgs = args;
            this.runtime.send(args);
        });
    }
    get state() {
        return this._state;
    }
    set state(state) {
        if (state === RuntimeWorkerState.BUSY) {
            this.socketReady = types_1.EmulatorLog.waitForLog(this.runtime.events, "SYSTEM", "runtime-status", (el) => {
                return el.data.state === "ready";
            });
        }
        if (state === RuntimeWorkerState.IDLE) {
            for (const l of this.logListeners) {
                this.runtime.events.removeListener("log", l);
            }
            this.logListeners = [];
            this.socketReady = undefined;
        }
        if (state === RuntimeWorkerState.FINISHED) {
            this.runtime.events.removeAllListeners();
        }
        this.log(state);
        this._state = state;
        this.stateEvents.emit(this._state);
    }
    onLogs(listener, forever = false) {
        if (!forever) {
            this.logListeners.push(listener);
        }
        this.runtime.events.on("log", listener);
    }
    waitForDone() {
        if (this.state === RuntimeWorkerState.IDLE || this.state === RuntimeWorkerState.FINISHED) {
            return Promise.resolve();
        }
        return new Promise((res) => {
            const listener = () => {
                this.stateEvents.removeListener(RuntimeWorkerState.IDLE, listener);
                this.stateEvents.removeListener(RuntimeWorkerState.FINISHED, listener);
                res();
            };
            this.stateEvents.once(RuntimeWorkerState.IDLE, listener);
            this.stateEvents.once(RuntimeWorkerState.FINISHED, listener);
        });
    }
    waitForSocketReady() {
        return (this.socketReady ||
            Promise.reject(new Error("Cannot call waitForSocketReady() if runtime is not BUSY")));
    }
    log(msg) {
        emulatorLogger_1.EmulatorLogger.log("DEBUG", `[worker-${this.key}-${this.id}]: ${msg}`);
    }
}
exports.RuntimeWorker = RuntimeWorker;
class RuntimeWorkerPool {
    constructor(mode = types_1.FunctionsExecutionMode.AUTO) {
        this.mode = mode;
        this.workers = new Map();
    }
    getKey(triggerId) {
        if (this.mode === types_1.FunctionsExecutionMode.SEQUENTIAL) {
            return "~shared~";
        }
        else {
            return triggerId || "~diagnostic~";
        }
    }
    refresh() {
        for (const arr of this.workers.values()) {
            arr.forEach((w) => {
                if (w.state === RuntimeWorkerState.IDLE) {
                    this.log(`Shutting down IDLE worker (${w.key})`);
                    w.state = RuntimeWorkerState.FINISHING;
                    w.runtime.shutdown();
                }
                else if (w.state === RuntimeWorkerState.BUSY) {
                    this.log(`Marking BUSY worker to finish (${w.key})`);
                    w.state = RuntimeWorkerState.FINISHING;
                }
            });
        }
    }
    exit() {
        for (const arr of this.workers.values()) {
            arr.forEach((w) => {
                if (w.state === RuntimeWorkerState.IDLE) {
                    w.runtime.shutdown();
                }
                else {
                    w.runtime.kill();
                }
            });
        }
    }
    readyForWork(triggerId) {
        const idleWorker = this.getIdleWorker(triggerId);
        return !!idleWorker;
    }
    submitWork(triggerId, frb, opts) {
        this.log(`submitWork(triggerId=${triggerId})`);
        const worker = this.getIdleWorker(triggerId);
        if (!worker) {
            throw new error_1.FirebaseError("Internal Error: can't call submitWork without checking for idle workers");
        }
        worker.execute(frb, opts);
        return worker;
    }
    getIdleWorker(triggerId) {
        this.cleanUpWorkers();
        const triggerWorkers = this.getTriggerWorkers(triggerId);
        if (!triggerWorkers.length) {
            this.setTriggerWorkers(triggerId, []);
            return;
        }
        for (const worker of triggerWorkers) {
            if (worker.state === RuntimeWorkerState.IDLE) {
                return worker;
            }
        }
        return;
    }
    addWorker(triggerId, runtime) {
        const worker = new RuntimeWorker(this.getKey(triggerId), runtime);
        this.log(`addWorker(${worker.key})`);
        const keyWorkers = this.getTriggerWorkers(triggerId);
        keyWorkers.push(worker);
        this.setTriggerWorkers(triggerId, keyWorkers);
        worker.onLogs((log) => {
            emulatorLogger_1.EmulatorLogger.handleRuntimeLog(log);
        }, true);
        this.log(`Adding worker with key ${worker.key}, total=${keyWorkers.length}`);
        return worker;
    }
    getTriggerWorkers(triggerId) {
        return this.workers.get(this.getKey(triggerId)) || [];
    }
    setTriggerWorkers(triggerId, workers) {
        this.workers.set(this.getKey(triggerId), workers);
    }
    cleanUpWorkers() {
        for (const [key, keyWorkers] of this.workers.entries()) {
            const notDoneWorkers = keyWorkers.filter((worker) => {
                return worker.state !== RuntimeWorkerState.FINISHED;
            });
            if (notDoneWorkers.length !== keyWorkers.length) {
                this.log(`Cleaned up workers for ${key}: ${keyWorkers.length} --> ${notDoneWorkers.length}`);
            }
            this.setTriggerWorkers(key, notDoneWorkers);
        }
    }
    log(msg) {
        emulatorLogger_1.EmulatorLogger.log("DEBUG", `[worker-pool] ${msg}`);
    }
}
exports.RuntimeWorkerPool = RuntimeWorkerPool;
