"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const API = require("./indexes-api");
const util = require("./util");
const QUERY_SCOPE_SEQUENCE = [
    API.QueryScope.COLLECTION_GROUP,
    API.QueryScope.COLLECTION,
    undefined,
];
const ORDER_SEQUENCE = [API.Order.ASCENDING, API.Order.DESCENDING, undefined];
const ARRAY_CONFIG_SEQUENCE = [API.ArrayConfig.CONTAINS, undefined];
function compareSpecIndex(a, b) {
    if (a.collectionGroup !== b.collectionGroup) {
        return a.collectionGroup.localeCompare(b.collectionGroup);
    }
    if (a.queryScope !== b.queryScope) {
        return compareQueryScope(a.queryScope, b.queryScope);
    }
    return compareArrays(a.fields, b.fields, compareIndexField);
}
exports.compareSpecIndex = compareSpecIndex;
function compareApiIndex(a, b) {
    if (a.name && b.name) {
        const aName = util.parseIndexName(a.name);
        const bName = util.parseIndexName(b.name);
        if (aName.collectionGroupId !== bName.collectionGroupId) {
            return aName.collectionGroupId.localeCompare(bName.collectionGroupId);
        }
    }
    if (a.queryScope !== b.queryScope) {
        return compareQueryScope(a.queryScope, b.queryScope);
    }
    return compareArrays(a.fields, b.fields, compareIndexField);
}
exports.compareApiIndex = compareApiIndex;
function compareApiField(a, b) {
    const aName = util.parseFieldName(a.name);
    const bName = util.parseFieldName(b.name);
    if (aName.collectionGroupId !== bName.collectionGroupId) {
        return aName.collectionGroupId.localeCompare(bName.collectionGroupId);
    }
    if (aName.fieldPath !== bName.fieldPath) {
        return aName.fieldPath.localeCompare(bName.fieldPath);
    }
    return compareArraysSorted(a.indexConfig.indexes || [], b.indexConfig.indexes || [], compareApiIndex);
}
exports.compareApiField = compareApiField;
function compareFieldOverride(a, b) {
    if (a.collectionGroup !== b.collectionGroup) {
        return a.collectionGroup.localeCompare(b.collectionGroup);
    }
    if (a.fieldPath !== b.fieldPath) {
        return a.fieldPath.localeCompare(b.fieldPath);
    }
    return compareArraysSorted(a.indexes, b.indexes, compareFieldIndex);
}
exports.compareFieldOverride = compareFieldOverride;
function compareIndexField(a, b) {
    if (a.fieldPath !== b.fieldPath) {
        return a.fieldPath.localeCompare(b.fieldPath);
    }
    if (a.order !== b.order) {
        return compareOrder(a.order, b.order);
    }
    if (a.arrayConfig !== b.arrayConfig) {
        return compareArrayConfig(a.arrayConfig, b.arrayConfig);
    }
    return 0;
}
function compareFieldIndex(a, b) {
    if (a.queryScope !== b.queryScope) {
        return compareQueryScope(a.queryScope, b.queryScope);
    }
    if (a.order !== b.order) {
        return compareOrder(a.order, b.order);
    }
    if (a.arrayConfig !== b.arrayConfig) {
        return compareArrayConfig(a.arrayConfig, b.arrayConfig);
    }
    return 0;
}
function compareQueryScope(a, b) {
    return QUERY_SCOPE_SEQUENCE.indexOf(a) - QUERY_SCOPE_SEQUENCE.indexOf(b);
}
function compareOrder(a, b) {
    return ORDER_SEQUENCE.indexOf(a) - ORDER_SEQUENCE.indexOf(b);
}
function compareArrayConfig(a, b) {
    return ARRAY_CONFIG_SEQUENCE.indexOf(a) - ARRAY_CONFIG_SEQUENCE.indexOf(b);
}
function compareArrays(a, b, fn) {
    const minFields = Math.min(a.length, b.length);
    for (let i = 0; i < minFields; i++) {
        const cmp = fn(a[i], b[i]);
        if (cmp !== 0) {
            return cmp;
        }
    }
    return a.length - b.length;
}
function compareArraysSorted(a, b, fn) {
    const aSorted = a.sort(fn);
    const bSorted = b.sort(fn);
    return compareArrays(aSorted, bSorted, fn);
}
