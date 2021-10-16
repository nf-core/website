"use strict";
var _ = require("lodash");
var clc = require("cli-color");
var marked = require("marked");
var { FirebaseError } = require("./error");
module.exports = function (options, allowNull = false) {
    if (!options.project && !allowNull) {
        var aliases = _.get(options, "rc.projects", {});
        var aliasCount = _.size(aliases);
        if (aliasCount === 0) {
            throw new FirebaseError("No currently active project.\n" +
                "To run this command, you need to specify a project. You have two options:\n" +
                "- Run this command with " +
                clc.bold("--project <alias_or_project_id>") +
                ".\n" +
                "- Set an active project by running " +
                clc.bold("firebase use --add") +
                ", then rerun this command.\n" +
                "To list all the Firebase projects to which you have access, run " +
                clc.bold("firebase projects:list") +
                ".\n" +
                marked("To learn about active projects for the CLI, visit https://firebase.google.com/docs/cli#project_aliases"), {
                exit: 1,
            });
        }
        else {
            var aliasList = _.map(aliases, function (projectId, aname) {
                return "  " + aname + " (" + projectId + ")";
            }).join("\n");
            throw new FirebaseError("No project active, but project aliases are available.\n\nRun " +
                clc.bold("firebase use <alias>") +
                " with one of these options:\n\n" +
                aliasList);
        }
    }
    return options.project;
};
