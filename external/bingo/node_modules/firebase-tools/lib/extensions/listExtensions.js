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
const _ = require("lodash");
const clc = require("cli-color");
const Table = require("cli-table");
const extensionsApi_1 = require("./extensionsApi");
const extensionsHelper_1 = require("./extensionsHelper");
const utils = require("../utils");
const logger = require("../logger");
function listExtensions(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        const instances = yield extensionsApi_1.listInstances(projectId);
        if (instances.length < 1) {
            utils.logLabeledBullet(extensionsHelper_1.logPrefix, `there are no extensions installed on project ${clc.bold(projectId)}.`);
            return { instances: [] };
        }
        const table = new Table({
            head: ["Extension Instance ID", "State", "Extension Version", "Create Time", "Update Time"],
            style: { head: ["yellow"] },
        });
        const sorted = _.sortBy(instances, "createTime", "asc").reverse();
        sorted.forEach((instance) => {
            table.push([
                _.last(instance.name.split("/")),
                instance.state,
                _.get(instance, "config.source.spec.version", ""),
                instance.createTime,
                instance.updateTime,
            ]);
        });
        utils.logLabeledBullet(extensionsHelper_1.logPrefix, `list of extensions installed in ${clc.bold(projectId)}:`);
        logger.info(table.toString());
        return { instances: sorted };
    });
}
exports.listExtensions = listExtensions;
