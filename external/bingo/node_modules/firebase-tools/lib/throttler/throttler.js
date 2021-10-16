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
const logger = require("../logger");
const retries_exhausted_error_1 = require("./errors/retries-exhausted-error");
const timeout_error_1 = require("./errors/timeout-error");
function backoff(retryNumber, delay) {
    return new Promise((resolve) => {
        setTimeout(resolve, delay * Math.pow(2, retryNumber));
    });
}
function DEFAULT_HANDLER(task) {
    return task();
}
class Throttler {
    constructor(options) {
        this.name = "";
        this.concurrency = 200;
        this.handler = DEFAULT_HANDLER;
        this.active = 0;
        this.complete = 0;
        this.success = 0;
        this.errored = 0;
        this.retried = 0;
        this.total = 0;
        this.taskDataMap = new Map();
        this.waits = [];
        this.min = 9999999999;
        this.max = 0;
        this.avg = 0;
        this.retries = 0;
        this.backoff = 200;
        this.closed = false;
        this.finished = false;
        this.startTime = 0;
        if (options.name) {
            this.name = options.name;
        }
        if (options.handler) {
            this.handler = options.handler;
        }
        if (typeof options.concurrency === "number") {
            this.concurrency = options.concurrency;
        }
        if (typeof options.retries === "number") {
            this.retries = options.retries;
        }
        if (typeof options.backoff === "number") {
            this.backoff = options.backoff;
        }
        if (typeof options.backoff === "number") {
            this.backoff = options.backoff;
        }
    }
    wait() {
        const p = new Promise((resolve, reject) => {
            this.waits.push({ resolve, reject });
        });
        return p;
    }
    add(task, timeoutMillis) {
        this.addHelper(task, timeoutMillis);
    }
    run(task, timeoutMillis) {
        return new Promise((resolve, reject) => {
            this.addHelper(task, timeoutMillis, { resolve, reject });
        });
    }
    close() {
        this.closed = true;
        return this.finishIfIdle();
    }
    process() {
        if (this.finishIfIdle() || this.active >= this.concurrency || !this.hasWaitingTask()) {
            return;
        }
        this.active++;
        this.handle(this.nextWaitingTaskIndex());
    }
    handle(cursorIndex) {
        return __awaiter(this, void 0, void 0, function* () {
            const taskData = this.taskDataMap.get(cursorIndex);
            if (!taskData) {
                throw new Error(`taskData.get(${cursorIndex}) does not exist`);
            }
            const promises = [this.executeTask(cursorIndex)];
            if (taskData.timeoutMillis) {
                promises.push(this.initializeTimeout(cursorIndex));
            }
            let result;
            try {
                result = yield Promise.race(promises);
            }
            catch (err) {
                this.errored++;
                this.complete++;
                this.active--;
                this.onTaskFailed(err, cursorIndex);
                return;
            }
            this.success++;
            this.complete++;
            this.active--;
            this.onTaskFulfilled(result, cursorIndex);
        });
    }
    stats() {
        return {
            max: this.max,
            min: this.min,
            avg: this.avg,
            active: this.active,
            complete: this.complete,
            success: this.success,
            errored: this.errored,
            retried: this.retried,
            total: this.total,
            elapsed: Date.now() - this.startTime,
        };
    }
    taskName(cursorIndex) {
        const taskData = this.taskDataMap.get(cursorIndex);
        if (!taskData) {
            return "finished task";
        }
        return typeof taskData.task === "string" ? taskData.task : `index ${cursorIndex}`;
    }
    addHelper(task, timeoutMillis, wait) {
        if (this.closed) {
            throw new Error("Cannot add a task to a closed throttler.");
        }
        if (!this.startTime) {
            this.startTime = Date.now();
        }
        this.taskDataMap.set(this.total, {
            task,
            wait,
            timeoutMillis,
            retryCount: 0,
            isTimedOut: false,
        });
        this.total++;
        this.process();
    }
    finishIfIdle() {
        if (this.closed && !this.hasWaitingTask() && this.active === 0) {
            this.finish();
            return true;
        }
        return false;
    }
    finish(err) {
        this.waits.forEach((p) => {
            if (err) {
                return p.reject(err);
            }
            this.finished = true;
            return p.resolve();
        });
    }
    initializeTimeout(cursorIndex) {
        const taskData = this.taskDataMap.get(cursorIndex);
        const timeoutMillis = taskData.timeoutMillis;
        const timeoutPromise = new Promise((_, reject) => {
            taskData.timeoutId = setTimeout(() => {
                taskData.isTimedOut = true;
                reject(new timeout_error_1.default(this.taskName(cursorIndex), timeoutMillis));
            }, timeoutMillis);
        });
        return timeoutPromise;
    }
    executeTask(cursorIndex) {
        return __awaiter(this, void 0, void 0, function* () {
            const taskData = this.taskDataMap.get(cursorIndex);
            const t0 = Date.now();
            let result;
            try {
                result = yield this.handler(taskData.task);
            }
            catch (err) {
                if (taskData.retryCount === this.retries) {
                    throw new retries_exhausted_error_1.default(this.taskName(cursorIndex), this.retries, err);
                }
                yield backoff(taskData.retryCount + 1, this.backoff);
                if (taskData.isTimedOut) {
                    throw new timeout_error_1.default(this.taskName(cursorIndex), taskData.timeoutMillis);
                }
                this.retried++;
                taskData.retryCount++;
                logger.debug(`[${this.name}] Retrying task`, this.taskName(cursorIndex));
                return this.executeTask(cursorIndex);
            }
            if (taskData.isTimedOut) {
                throw new timeout_error_1.default(this.taskName(cursorIndex), taskData.timeoutMillis);
            }
            const dt = Date.now() - t0;
            this.min = Math.min(dt, this.min);
            this.max = Math.max(dt, this.max);
            this.avg = (this.avg * this.complete + dt) / (this.complete + 1);
            return result;
        });
    }
    onTaskFulfilled(result, cursorIndex) {
        const taskData = this.taskDataMap.get(cursorIndex);
        if (taskData.wait) {
            taskData.wait.resolve(result);
        }
        this.cleanupTask(cursorIndex);
        this.process();
    }
    onTaskFailed(error, cursorIndex) {
        const taskData = this.taskDataMap.get(cursorIndex);
        logger.debug(error);
        if (taskData.wait) {
            taskData.wait.reject(error);
        }
        this.cleanupTask(cursorIndex);
        this.finish(error);
    }
    cleanupTask(cursorIndex) {
        const { timeoutId } = this.taskDataMap.get(cursorIndex);
        if (timeoutId) {
            clearTimeout(timeoutId);
        }
        this.taskDataMap.delete(cursorIndex);
    }
}
exports.Throttler = Throttler;
