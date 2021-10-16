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
const emulatorLogger_1 = require("./emulatorLogger");
const types_1 = require("./types");
class WorkQueue {
    constructor(mode = types_1.FunctionsExecutionMode.AUTO) {
        this.mode = mode;
        this.queue = [];
        this.workRunningCount = 0;
        this.notifyQueue = () => { };
        this.stopped = true;
    }
    submit(entry) {
        this.queue.push(entry);
        this.notifyQueue();
        this.logState();
    }
    start() {
        return __awaiter(this, void 0, void 0, function* () {
            if (!this.stopped) {
                return;
            }
            this.stopped = false;
            while (!this.stopped) {
                if (!this.queue.length) {
                    yield new Promise((res) => {
                        this.notifyQueue = res;
                    });
                }
                const workPromise = this.runNext();
                if (this.mode === types_1.FunctionsExecutionMode.SEQUENTIAL) {
                    yield workPromise;
                }
            }
        });
    }
    stop() {
        this.stopped = true;
    }
    runNext() {
        return __awaiter(this, void 0, void 0, function* () {
            const next = this.queue.shift();
            if (next) {
                this.workRunningCount++;
                this.logState();
                try {
                    yield next();
                }
                catch (e) {
                    emulatorLogger_1.EmulatorLogger.log("DEBUG", e);
                }
                finally {
                    this.workRunningCount--;
                    this.logState();
                }
            }
        });
    }
    logState() {
        emulatorLogger_1.EmulatorLogger.logLabeled("DEBUG", "work-queue", JSON.stringify({
            queueLength: this.queue.length,
            workRunningCount: this.workRunningCount,
        }));
    }
}
exports.WorkQueue = WorkQueue;
