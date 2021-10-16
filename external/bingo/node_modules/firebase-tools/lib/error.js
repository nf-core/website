"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const lodash_1 = require("lodash");
const DEFAULT_CHILDREN = [];
const DEFAULT_EXIT = 1;
const DEFAULT_STATUS = 500;
class FirebaseError extends Error {
    constructor(message, options = {}) {
        super();
        this.name = "FirebaseError";
        this.children = lodash_1.defaultTo(options.children, DEFAULT_CHILDREN);
        this.context = options.context;
        this.exit = lodash_1.defaultTo(options.exit, DEFAULT_EXIT);
        this.message = message;
        this.original = options.original;
        this.status = lodash_1.defaultTo(options.status, DEFAULT_STATUS);
    }
}
exports.FirebaseError = FirebaseError;
