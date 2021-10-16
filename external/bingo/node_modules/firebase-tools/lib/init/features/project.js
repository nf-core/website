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
const _ = require("lodash");
const error_1 = require("../../error");
const projects_1 = require("../../management/projects");
const logger = require("../../logger");
const prompt_1 = require("../../prompt");
const utils = require("../../utils");
const OPTION_NO_PROJECT = "Don't set up a default project";
const OPTION_USE_PROJECT = "Use an existing project";
const OPTION_NEW_PROJECT = "Create a new project";
const OPTION_ADD_FIREBASE = "Add Firebase to an existing Google Cloud Platform project";
function toProjectInfo(projectMetaData) {
    const { projectId, displayName, resources } = projectMetaData;
    return {
        id: projectId,
        label: `${projectId}` + (displayName ? ` (${displayName})` : ""),
        instance: _.get(resources, "realtimeDatabaseInstance"),
        location: _.get(resources, "locationId"),
    };
}
function promptAndCreateNewProject() {
    return __awaiter(this, void 0, void 0, function* () {
        utils.logBullet("If you want to create a project in a Google Cloud organization or folder, please use " +
            `"firebase projects:create" instead, and return to this command when you've created the project.`);
        const promptAnswer = {};
        yield prompt_1.prompt(promptAnswer, projects_1.PROJECTS_CREATE_QUESTIONS);
        if (!promptAnswer.projectId) {
            throw new error_1.FirebaseError("Project ID cannot be empty");
        }
        return yield projects_1.createFirebaseProjectAndLog(promptAnswer.projectId, {
            displayName: promptAnswer.displayName,
        });
    });
}
function promptAndAddFirebaseToCloudProject() {
    return __awaiter(this, void 0, void 0, function* () {
        const projectId = yield projects_1.promptAvailableProjectId();
        if (!projectId) {
            throw new error_1.FirebaseError("Project ID cannot be empty");
        }
        return yield projects_1.addFirebaseToCloudProjectAndLog(projectId);
    });
}
function projectChoicePrompt(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const choices = [
            { name: OPTION_USE_PROJECT, value: OPTION_USE_PROJECT },
            { name: OPTION_NEW_PROJECT, value: OPTION_NEW_PROJECT },
            { name: OPTION_ADD_FIREBASE, value: OPTION_ADD_FIREBASE },
            { name: OPTION_NO_PROJECT, value: OPTION_NO_PROJECT },
        ];
        const projectSetupOption = yield prompt_1.promptOnce({
            type: "list",
            name: "id",
            message: "Please select an option:",
            choices,
        });
        switch (projectSetupOption) {
            case OPTION_USE_PROJECT:
                return projects_1.getOrPromptProject(options);
            case OPTION_NEW_PROJECT:
                return promptAndCreateNewProject();
            case OPTION_ADD_FIREBASE:
                return promptAndAddFirebaseToCloudProject();
            default:
                return;
        }
    });
}
function doSetup(setup, config, options) {
    return __awaiter(this, void 0, void 0, function* () {
        setup.project = {};
        logger.info();
        logger.info(`First, let's associate this project directory with a Firebase project.`);
        logger.info(`You can create multiple project aliases by running ${clc.bold("firebase use --add")}, `);
        logger.info(`but for now we'll just set up a default project.`);
        logger.info();
        const projectFromRcFile = _.get(setup.rcfile, "projects.default");
        if (projectFromRcFile) {
            utils.logBullet(`.firebaserc already has a default project, using ${projectFromRcFile}.`);
            const rcProject = yield projects_1.getFirebaseProject(projectFromRcFile);
            setup.projectId = rcProject.projectId;
            setup.projectLocation = _.get(rcProject, "resources.locationId");
            return;
        }
        let projectMetaData;
        if (options.project) {
            projectMetaData = yield projects_1.getFirebaseProject(options.project);
        }
        else {
            projectMetaData = yield projectChoicePrompt(options);
            if (!projectMetaData) {
                return;
            }
        }
        const projectInfo = toProjectInfo(projectMetaData);
        utils.logBullet(`Using project ${projectInfo.label}`);
        _.set(setup.rcfile, "projects.default", projectInfo.id);
        setup.projectId = projectInfo.id;
        setup.instance = projectInfo.instance;
        setup.projectLocation = projectInfo.location;
        utils.makeActiveProject(config.projectDir, projectInfo.id);
    });
}
exports.doSetup = doSetup;
