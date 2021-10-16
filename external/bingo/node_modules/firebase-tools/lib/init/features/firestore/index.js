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
const ensureCloudResourceLocation_1 = require("../../../ensureCloudResourceLocation");
const requirePermissions_1 = require("../../../requirePermissions");
const rules = require("./rules");
const indexes = require("./indexes");
function doSetup(setup, config) {
    return __awaiter(this, void 0, void 0, function* () {
        setup.config.firestore = {};
        ensureCloudResourceLocation_1.ensureLocationSet(setup.projectLocation, "Cloud Firestore");
        yield requirePermissions_1.requirePermissions({ project: setup.projectId });
        yield rules.initRules(setup, config);
        yield indexes.initIndexes(setup, config);
    });
}
exports.doSetup = doSetup;
