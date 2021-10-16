"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var Mode;
(function (Mode) {
    Mode["ASCENDING"] = "ASCENDING";
    Mode["DESCENDING"] = "DESCENDING";
    Mode["ARRAY_CONTAINS"] = "ARRAY_CONTAINS";
})(Mode = exports.Mode || (exports.Mode = {}));
var QueryScope;
(function (QueryScope) {
    QueryScope["COLLECTION"] = "COLLECTION";
    QueryScope["COLLECTION_GROUP"] = "COLLECTION_GROUP";
})(QueryScope = exports.QueryScope || (exports.QueryScope = {}));
var Order;
(function (Order) {
    Order["ASCENDING"] = "ASCENDING";
    Order["DESCENDING"] = "DESCENDING";
})(Order = exports.Order || (exports.Order = {}));
var ArrayConfig;
(function (ArrayConfig) {
    ArrayConfig["CONTAINS"] = "CONTAINS";
})(ArrayConfig = exports.ArrayConfig || (exports.ArrayConfig = {}));
var State;
(function (State) {
    State["CREATING"] = "CREATING";
    State["READY"] = "READY";
    State["NEEDS_REPAIR"] = "NEEDS_REPAIR";
})(State = exports.State || (exports.State = {}));
