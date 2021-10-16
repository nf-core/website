"use strict";
var _ = require("lodash");
var api = require("./api");
module.exports = function (id) {
    return api.getProjects().then(function (projects) {
        if (_.includes(_.keys(projects), id)) {
            return id;
        }
        for (var projectId in projects) {
            if (projects.hasOwnProperty(projectId)) {
                var instance = _.get(projects, [projectId, "instances", "database", "0"]);
                if (id === instance) {
                    return projectId;
                }
            }
        }
        return null;
    });
};
