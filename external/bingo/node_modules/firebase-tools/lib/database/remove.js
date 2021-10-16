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
const pathLib = require("path");
const removeRemote_1 = require("./removeRemote");
const listRemote_1 = require("./listRemote");
const stack_1 = require("../throttler/stack");
function chunkList(ls, chunkSize) {
    const chunks = [];
    for (let i = 0; i < ls.length; i += chunkSize) {
        chunks.push(ls.slice(i, i + chunkSize));
    }
    return chunks;
}
const INITIAL_DELETE_BATCH_SIZE = 25;
const INITIAL_LIST_NUM_SUB_PATH = 100;
const MAX_LIST_NUM_SUB_PATH = 204800;
class DatabaseRemove {
    constructor(instance, path) {
        this.path = path;
        this.remote = new removeRemote_1.RTDBRemoveRemote(instance);
        this.deleteJobStack = new stack_1.Stack({
            name: "delete stack",
            concurrency: 1,
            retries: 3,
        });
        this.listRemote = new listRemote_1.RTDBListRemote(instance);
        this.listStack = new stack_1.Stack({
            name: "list stack",
            concurrency: 1,
            retries: 3,
        });
    }
    execute() {
        return __awaiter(this, void 0, void 0, function* () {
            yield this.deletePath(this.path);
        });
    }
    deletePath(path) {
        return __awaiter(this, void 0, void 0, function* () {
            if (yield this.deleteJobStack.run(() => this.remote.deletePath(path))) {
                return Promise.resolve(true);
            }
            let listNumSubPath = INITIAL_LIST_NUM_SUB_PATH;
            let batchSizeLow = 1;
            let batchSizeHigh = MAX_LIST_NUM_SUB_PATH + 1;
            let batchSize = INITIAL_DELETE_BATCH_SIZE;
            while (true) {
                const subPathList = yield this.listStack.run(() => this.listRemote.listPath(path, listNumSubPath));
                if (subPathList.length === 0) {
                    return Promise.resolve(false);
                }
                const chunks = chunkList(subPathList, batchSize);
                let nSmallChunks = 0;
                for (const chunk of chunks) {
                    if (yield this.deleteSubPath(path, chunk)) {
                        nSmallChunks += 1;
                    }
                }
                if (nSmallChunks > chunks.length / 2) {
                    batchSizeLow = batchSize;
                    batchSize = Math.floor(Math.min(batchSize * 2, (batchSizeHigh + batchSize) / 2));
                }
                else {
                    batchSizeHigh = batchSize;
                    batchSize = Math.floor((batchSizeLow + batchSize) / 2);
                }
                if (listNumSubPath * 2 <= MAX_LIST_NUM_SUB_PATH) {
                    listNumSubPath = listNumSubPath * 2;
                }
                else {
                    listNumSubPath = Math.floor(MAX_LIST_NUM_SUB_PATH / batchSize) * batchSize;
                }
            }
        });
    }
    deleteSubPath(path, subPaths) {
        return __awaiter(this, void 0, void 0, function* () {
            if (subPaths.length === 0) {
                throw new Error("deleteSubPath is called with empty subPaths list");
            }
            if (subPaths.length === 1) {
                return this.deletePath(pathLib.join(path, subPaths[0]));
            }
            if (yield this.deleteJobStack.run(() => this.remote.deleteSubPath(path, subPaths))) {
                return Promise.resolve(true);
            }
            const mid = Math.floor(subPaths.length / 2);
            yield this.deleteSubPath(path, subPaths.slice(0, mid));
            yield this.deleteSubPath(path, subPaths.slice(mid));
            return Promise.resolve(false);
        });
    }
}
exports.default = DatabaseRemove;
