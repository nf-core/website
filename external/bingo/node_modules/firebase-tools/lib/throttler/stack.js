"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const throttler_1 = require("./throttler");
class Stack extends throttler_1.Throttler {
    constructor(options) {
        super(options);
        this.lastTotal = 0;
        this.stack = [];
        this.name = this.name || "stack";
    }
    hasWaitingTask() {
        return this.lastTotal !== this.total || this.stack.length > 0;
    }
    nextWaitingTaskIndex() {
        while (this.lastTotal < this.total) {
            this.stack.push(this.lastTotal);
            this.lastTotal++;
        }
        const next = this.stack.pop();
        if (next === undefined) {
            throw new Error("There is no more task in stack");
        }
        return next;
    }
}
exports.Stack = Stack;
exports.default = Stack;
