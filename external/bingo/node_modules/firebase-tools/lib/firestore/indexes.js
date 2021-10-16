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
const clc = require("cli-color");
const api = require("../api");
const logger = require("../logger");
const utils = require("../utils");
const validator = require("./validator");
const API = require("./indexes-api");
const sort = require("./indexes-sort");
const util = require("./util");
class FirestoreIndexes {
    deploy(project, indexes, fieldOverrides) {
        return __awaiter(this, void 0, void 0, function* () {
            const spec = this.upgradeOldSpec({
                indexes,
                fieldOverrides,
            });
            this.validateSpec(spec);
            const indexesToDeploy = spec.indexes;
            const fieldOverridesToDeploy = spec.fieldOverrides;
            const existingIndexes = yield this.listIndexes(project);
            const existingFieldOverrides = yield this.listFieldOverrides(project);
            if (existingIndexes.length > indexesToDeploy.length) {
                utils.logBullet(clc.bold.cyan("firestore:") +
                    " there are some indexes defined in your project that are not present in your " +
                    "firestore indexes file. Run firebase firestore:indexes and save the result to correct the discrepancy.");
            }
            for (const index of indexesToDeploy) {
                const exists = existingIndexes.some((x) => this.indexMatchesSpec(x, index));
                if (exists) {
                    logger.debug(`Skipping existing index: ${JSON.stringify(index)}`);
                }
                else {
                    logger.debug(`Creating new index: ${JSON.stringify(index)}`);
                    yield this.createIndex(project, index);
                }
            }
            if (existingFieldOverrides.length > fieldOverridesToDeploy.length) {
                utils.logBullet(clc.bold.cyan("firestore:") +
                    " there are some field overrides defined in your project that are not present in your " +
                    "firestore indexes file. Run firebase firestore:indexes and save the result to correct the discrepancy.");
            }
            for (const field of fieldOverridesToDeploy) {
                const exists = existingFieldOverrides.some((x) => this.fieldMatchesSpec(x, field));
                if (exists) {
                    logger.debug(`Skipping existing field override: ${JSON.stringify(field)}`);
                }
                else {
                    logger.debug(`Updating field override: ${JSON.stringify(field)}`);
                    yield this.patchField(project, field);
                }
            }
        });
    }
    listIndexes(project) {
        return __awaiter(this, void 0, void 0, function* () {
            const url = `projects/${project}/databases/(default)/collectionGroups/-/indexes`;
            const res = yield api.request("GET", `/v1beta2/${url}`, {
                auth: true,
                origin: api.firestoreOrigin,
            });
            const indexes = res.body.indexes;
            if (!indexes) {
                return [];
            }
            return indexes.map((index) => {
                const fields = index.fields.filter((field) => {
                    return field.fieldPath !== "__name__";
                });
                return {
                    name: index.name,
                    state: index.state,
                    queryScope: index.queryScope,
                    fields,
                };
            });
        });
    }
    listFieldOverrides(project) {
        return __awaiter(this, void 0, void 0, function* () {
            const parent = `projects/${project}/databases/(default)/collectionGroups/-`;
            const url = `${parent}/fields?filter=indexConfig.usesAncestorConfig=false`;
            const res = yield api.request("GET", `/v1beta2/${url}`, {
                auth: true,
                origin: api.firestoreOrigin,
            });
            const fields = res.body.fields;
            if (!fields) {
                return [];
            }
            return fields.filter((field) => {
                return field.name.indexOf("__default__") < 0;
            });
        });
    }
    makeIndexSpec(indexes, fields) {
        const indexesJson = indexes.map((index) => {
            return {
                collectionGroup: util.parseIndexName(index.name).collectionGroupId,
                queryScope: index.queryScope,
                fields: index.fields,
            };
        });
        if (!fields) {
            logger.debug("No field overrides specified, using [].");
            fields = [];
        }
        const fieldsJson = fields.map((field) => {
            const parsedName = util.parseFieldName(field.name);
            const fieldIndexes = field.indexConfig.indexes || [];
            return {
                collectionGroup: parsedName.collectionGroupId,
                fieldPath: parsedName.fieldPath,
                indexes: fieldIndexes.map((index) => {
                    const firstField = index.fields[0];
                    return {
                        order: firstField.order,
                        arrayConfig: firstField.arrayConfig,
                        queryScope: index.queryScope,
                    };
                }),
            };
        });
        const sortedIndexes = indexesJson.sort(sort.compareSpecIndex);
        const sortedFields = fieldsJson.sort(sort.compareFieldOverride);
        return {
            indexes: sortedIndexes,
            fieldOverrides: sortedFields,
        };
    }
    prettyPrintIndexes(indexes) {
        if (indexes.length === 0) {
            logger.info("None");
            return;
        }
        const sortedIndexes = indexes.sort(sort.compareApiIndex);
        sortedIndexes.forEach((index) => {
            logger.info(this.prettyIndexString(index));
        });
    }
    printFieldOverrides(fields) {
        if (fields.length === 0) {
            logger.info("None");
            return;
        }
        const sortedFields = fields.sort(sort.compareApiField);
        sortedFields.forEach((field) => {
            logger.info(this.prettyFieldString(field));
        });
    }
    validateSpec(spec) {
        validator.assertHas(spec, "indexes");
        spec.indexes.forEach((index) => {
            this.validateIndex(index);
        });
        if (spec.fieldOverrides) {
            spec.fieldOverrides.forEach((field) => {
                this.validateField(field);
            });
        }
    }
    validateIndex(index) {
        validator.assertHas(index, "collectionGroup");
        validator.assertHas(index, "queryScope");
        validator.assertEnum(index, "queryScope", Object.keys(API.QueryScope));
        validator.assertHas(index, "fields");
        index.fields.forEach((field) => {
            validator.assertHas(field, "fieldPath");
            validator.assertHasOneOf(field, ["order", "arrayConfig"]);
            if (field.order) {
                validator.assertEnum(field, "order", Object.keys(API.Order));
            }
            if (field.arrayConfig) {
                validator.assertEnum(field, "arrayConfig", Object.keys(API.ArrayConfig));
            }
        });
    }
    validateField(field) {
        validator.assertHas(field, "collectionGroup");
        validator.assertHas(field, "fieldPath");
        validator.assertHas(field, "indexes");
        field.indexes.forEach((index) => {
            validator.assertHasOneOf(index, ["arrayConfig", "order"]);
            if (index.arrayConfig) {
                validator.assertEnum(index, "arrayConfig", Object.keys(API.ArrayConfig));
            }
            if (index.order) {
                validator.assertEnum(index, "order", Object.keys(API.Order));
            }
            if (index.queryScope) {
                validator.assertEnum(index, "queryScope", Object.keys(API.QueryScope));
            }
        });
    }
    patchField(project, spec) {
        return __awaiter(this, void 0, void 0, function* () {
            const url = `projects/${project}/databases/(default)/collectionGroups/${spec.collectionGroup}/fields/${spec.fieldPath}`;
            const indexes = spec.indexes.map((index) => {
                return {
                    queryScope: index.queryScope,
                    fields: [
                        {
                            fieldPath: spec.fieldPath,
                            arrayConfig: index.arrayConfig,
                            order: index.order,
                        },
                    ],
                };
            });
            const data = {
                indexConfig: {
                    indexes,
                },
            };
            yield api.request("PATCH", `/v1beta2/${url}`, {
                auth: true,
                origin: api.firestoreOrigin,
                data,
            });
        });
    }
    createIndex(project, index) {
        const url = `projects/${project}/databases/(default)/collectionGroups/${index.collectionGroup}/indexes`;
        return api.request("POST", "/v1beta2/" + url, {
            auth: true,
            data: {
                fields: index.fields,
                queryScope: index.queryScope,
            },
            origin: api.firestoreOrigin,
        });
    }
    indexMatchesSpec(index, spec) {
        const collection = util.parseIndexName(index.name).collectionGroupId;
        if (collection !== spec.collectionGroup) {
            return false;
        }
        if (index.queryScope !== spec.queryScope) {
            return false;
        }
        if (index.fields.length !== spec.fields.length) {
            return false;
        }
        let i = 0;
        while (i < index.fields.length) {
            const iField = index.fields[i];
            const sField = spec.fields[i];
            if (iField.fieldPath !== sField.fieldPath) {
                return false;
            }
            if (iField.order !== sField.order) {
                return false;
            }
            if (iField.arrayConfig !== sField.arrayConfig) {
                return false;
            }
            i++;
        }
        return true;
    }
    fieldMatchesSpec(field, spec) {
        const parsedName = util.parseFieldName(field.name);
        if (parsedName.collectionGroupId !== spec.collectionGroup) {
            return false;
        }
        if (parsedName.fieldPath !== spec.fieldPath) {
            return false;
        }
        const fieldIndexes = field.indexConfig.indexes || [];
        if (fieldIndexes.length !== spec.indexes.length) {
            return false;
        }
        const fieldModes = fieldIndexes.map((index) => {
            const firstField = index.fields[0];
            return firstField.order || firstField.arrayConfig;
        });
        const specModes = spec.indexes.map((index) => {
            return index.order || index.arrayConfig;
        });
        for (const mode of fieldModes) {
            if (specModes.indexOf(mode) < 0) {
                return false;
            }
        }
        return true;
    }
    upgradeOldSpec(spec) {
        const result = {
            indexes: [],
            fieldOverrides: spec.fieldOverrides || [],
        };
        if (!(spec.indexes && spec.indexes.length > 0)) {
            return result;
        }
        if (spec.indexes[0].collectionId) {
            utils.logBullet(clc.bold.cyan("firestore:") +
                " your indexes indexes are specified in the v1beta1 API format. " +
                "Please upgrade to the new index API format by running " +
                clc.bold("firebase firestore:indexes") +
                " again and saving the result.");
        }
        result.indexes = spec.indexes.map((index) => {
            const i = {
                collectionGroup: index.collectionGroup || index.collectionId,
                queryScope: index.queryScope || API.QueryScope.COLLECTION,
                fields: [],
            };
            if (index.fields) {
                i.fields = index.fields.map((field) => {
                    const f = {
                        fieldPath: field.fieldPath,
                    };
                    if (field.order) {
                        f.order = field.order;
                    }
                    else if (field.arrayConfig) {
                        f.arrayConfig = field.arrayConfig;
                    }
                    else if (field.mode === API.Mode.ARRAY_CONTAINS) {
                        f.arrayConfig = API.ArrayConfig.CONTAINS;
                    }
                    else {
                        f.order = field.mode;
                    }
                    return f;
                });
            }
            return i;
        });
        return result;
    }
    prettyIndexString(index) {
        let result = "";
        if (index.state) {
            const stateMsg = `[${index.state}] `;
            if (index.state === API.State.READY) {
                result += clc.green(stateMsg);
            }
            else if (index.state === API.State.CREATING) {
                result += clc.yellow(stateMsg);
            }
            else {
                result += clc.red(stateMsg);
            }
        }
        const nameInfo = util.parseIndexName(index.name);
        result += clc.cyan(`(${nameInfo.collectionGroupId})`);
        result += " -- ";
        index.fields.forEach((field) => {
            if (field.fieldPath === "__name__") {
                return;
            }
            const orderOrArrayConfig = field.order ? field.order : field.arrayConfig;
            result += `(${field.fieldPath},${orderOrArrayConfig}) `;
        });
        return result;
    }
    prettyFieldString(field) {
        let result = "";
        const parsedName = util.parseFieldName(field.name);
        result +=
            "[" +
                clc.cyan(parsedName.collectionGroupId) +
                "." +
                clc.yellow(parsedName.fieldPath) +
                "] --";
        const fieldIndexes = field.indexConfig.indexes || [];
        if (fieldIndexes.length > 0) {
            fieldIndexes.forEach((index) => {
                const firstField = index.fields[0];
                const mode = firstField.order || firstField.arrayConfig;
                result += ` (${mode})`;
            });
        }
        else {
            result += " (no indexes)";
        }
        return result;
    }
}
exports.FirestoreIndexes = FirestoreIndexes;
