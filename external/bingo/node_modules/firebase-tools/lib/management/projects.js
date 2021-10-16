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
const api = require("../api");
const clc = require("cli-color");
const ora = require("ora");
const _ = require("lodash");
const logger = require("../logger");
const error_1 = require("../error");
const operation_poller_1 = require("../operation-poller");
const prompt_1 = require("../prompt");
const utils = require("../utils");
const TIMEOUT_MILLIS = 30000;
const MAXIMUM_PROMPT_LIST = 100;
const PROJECT_LIST_PAGE_SIZE = 1000;
const CREATE_PROJECT_API_REQUEST_TIMEOUT_MILLIS = 15000;
var ProjectParentResourceType;
(function (ProjectParentResourceType) {
    ProjectParentResourceType["ORGANIZATION"] = "organization";
    ProjectParentResourceType["FOLDER"] = "folder";
})(ProjectParentResourceType = exports.ProjectParentResourceType || (exports.ProjectParentResourceType = {}));
exports.PROJECTS_CREATE_QUESTIONS = [
    {
        type: "input",
        name: "projectId",
        default: "",
        message: "Please specify a unique project id " +
            `(${clc.yellow("warning")}: cannot be modified afterward) [6-30 characters]:\n`,
    },
    {
        type: "input",
        name: "displayName",
        default: "",
        message: "What would you like to call your project? (defaults to your project ID)",
    },
];
function createFirebaseProjectAndLog(projectId, options) {
    return __awaiter(this, void 0, void 0, function* () {
        const spinner = ora("Creating Google Cloud Platform project").start();
        try {
            yield createCloudProject(projectId, options);
            spinner.succeed();
        }
        catch (err) {
            spinner.fail();
            throw err;
        }
        return addFirebaseToCloudProjectAndLog(projectId);
    });
}
exports.createFirebaseProjectAndLog = createFirebaseProjectAndLog;
function addFirebaseToCloudProjectAndLog(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        let projectInfo;
        const spinner = ora("Adding Firebase resources to Google Cloud Platform project").start();
        try {
            projectInfo = yield addFirebaseToCloudProject(projectId);
        }
        catch (err) {
            spinner.fail();
            throw err;
        }
        spinner.succeed();
        logNewFirebaseProjectInfo(projectInfo);
        return projectInfo;
    });
}
exports.addFirebaseToCloudProjectAndLog = addFirebaseToCloudProjectAndLog;
function logNewFirebaseProjectInfo(projectInfo) {
    logger.info("");
    if (process.platform === "win32") {
        logger.info("=== Your Firebase project is ready! ===");
    }
    else {
        logger.info("🎉🎉🎉 Your Firebase project is ready! 🎉🎉🎉");
    }
    logger.info("");
    logger.info("Project information:");
    logger.info(`   - Project ID: ${clc.bold(projectInfo.projectId)}`);
    logger.info(`   - Project Name: ${clc.bold(projectInfo.displayName)}`);
    logger.info("");
    logger.info("Firebase console is available at");
    logger.info(`https://console.firebase.google.com/project/${clc.bold(projectInfo.projectId)}/overview`);
}
function getOrPromptProject(options) {
    return __awaiter(this, void 0, void 0, function* () {
        if (options.project) {
            return yield getFirebaseProject(options.project);
        }
        return selectProjectInteractively();
    });
}
exports.getOrPromptProject = getOrPromptProject;
function selectProjectInteractively(pageSize = MAXIMUM_PROMPT_LIST) {
    return __awaiter(this, void 0, void 0, function* () {
        const { projects, nextPageToken } = yield getFirebaseProjectPage(pageSize);
        if (projects.length === 0) {
            throw new error_1.FirebaseError("There are no Firebase projects associated with this account.");
        }
        if (nextPageToken) {
            return selectProjectByPrompting();
        }
        return selectProjectFromList(projects);
    });
}
function selectProjectByPrompting() {
    return __awaiter(this, void 0, void 0, function* () {
        const projectId = yield prompt_1.promptOnce({
            type: "input",
            message: "Please input the project ID you would like to use:",
        });
        return yield getFirebaseProject(projectId);
    });
}
function selectProjectFromList(projects = []) {
    return __awaiter(this, void 0, void 0, function* () {
        let choices = projects
            .filter((p) => !!p)
            .map((p) => {
            return {
                name: p.projectId + (p.displayName ? ` (${p.displayName})` : ""),
                value: p.projectId,
            };
        });
        choices = _.orderBy(choices, ["name"], ["asc"]);
        if (choices.length >= 25) {
            utils.logBullet(`Don't want to scroll through all your projects? If you know your project ID, ` +
                `you can initialize it directly using ${clc.bold("firebase init --project <project_id>")}.\n`);
        }
        const projectId = yield prompt_1.promptOnce({
            type: "list",
            name: "id",
            message: "Select a default Firebase project for this directory:",
            choices,
        });
        const project = projects.find((p) => p.projectId === projectId);
        if (!project) {
            throw new error_1.FirebaseError("Unexpected error. Project does not exist");
        }
        return project;
    });
}
function getProjectId(cloudProject) {
    const resourceName = cloudProject.project;
    return resourceName.substring(resourceName.lastIndexOf("/") + 1);
}
function promptAvailableProjectId() {
    return __awaiter(this, void 0, void 0, function* () {
        const { projects, nextPageToken } = yield getAvailableCloudProjectPage(MAXIMUM_PROMPT_LIST);
        if (projects.length === 0) {
            throw new error_1.FirebaseError("There are no available Google Cloud projects to add Firebase services.");
        }
        if (nextPageToken) {
            return yield prompt_1.promptOnce({
                type: "input",
                message: "Please input the ID of the Google Cloud Project you would like to add Firebase:",
            });
        }
        else {
            let choices = projects
                .filter((p) => !!p)
                .map((p) => {
                const projectId = getProjectId(p);
                return {
                    name: projectId + (p.displayName ? ` (${p.displayName})` : ""),
                    value: projectId,
                };
            });
            choices = _.orderBy(choices, ["name"], ["asc"]);
            return yield prompt_1.promptOnce({
                type: "list",
                name: "id",
                message: "Select the Google Cloud Platform project you would like to add Firebase:",
                choices,
            });
        }
    });
}
exports.promptAvailableProjectId = promptAvailableProjectId;
function createCloudProject(projectId, options) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const response = yield api.request("POST", "/v1/projects", {
                auth: true,
                origin: api.resourceManagerOrigin,
                timeout: CREATE_PROJECT_API_REQUEST_TIMEOUT_MILLIS,
                data: { projectId, name: options.displayName || projectId, parent: options.parentResource },
            });
            const projectInfo = yield operation_poller_1.pollOperation({
                pollerName: "Project Creation Poller",
                apiOrigin: api.resourceManagerOrigin,
                apiVersion: "v1",
                operationResourceName: response.body.name,
            });
            return projectInfo;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError("Failed to create Google Cloud project. See firebase-debug.log for more info.", { exit: 2, original: err });
        }
    });
}
exports.createCloudProject = createCloudProject;
function addFirebaseToCloudProject(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const response = yield api.request("POST", `/v1beta1/projects/${projectId}:addFirebase`, {
                auth: true,
                origin: api.firebaseApiOrigin,
                timeout: CREATE_PROJECT_API_REQUEST_TIMEOUT_MILLIS,
            });
            const projectInfo = yield operation_poller_1.pollOperation({
                pollerName: "Add Firebase Poller",
                apiOrigin: api.firebaseApiOrigin,
                apiVersion: "v1beta1",
                operationResourceName: response.body.name,
            });
            return projectInfo;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError("Failed to add Firebase to Google Cloud Platform project. See firebase-debug.log for more info.", { exit: 2, original: err });
        }
    });
}
exports.addFirebaseToCloudProject = addFirebaseToCloudProject;
function getProjectPage(apiResource, options) {
    return __awaiter(this, void 0, void 0, function* () {
        let apiResponse;
        const { responseKey, pageToken, pageSize } = options;
        const pageTokenQueryString = pageToken ? `&pageToken=${pageToken}` : "";
        apiResponse = yield api.request("GET", `${apiResource}?pageSize=${pageSize}${pageTokenQueryString}`, {
            auth: true,
            origin: api.firebaseApiOrigin,
            timeout: TIMEOUT_MILLIS,
        });
        return {
            projects: apiResponse.body[responseKey] || [],
            nextPageToken: apiResponse.body.nextPageToken,
        };
    });
}
function getFirebaseProjectPage(pageSize = PROJECT_LIST_PAGE_SIZE, pageToken) {
    return __awaiter(this, void 0, void 0, function* () {
        let projectPage;
        try {
            projectPage = yield getProjectPage("/v1beta1/projects", {
                responseKey: "results",
                pageSize,
                pageToken,
            });
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError("Failed to list Firebase projects. See firebase-debug.log for more info.", { exit: 2, original: err });
        }
        return projectPage;
    });
}
exports.getFirebaseProjectPage = getFirebaseProjectPage;
function getAvailableCloudProjectPage(pageSize = PROJECT_LIST_PAGE_SIZE, pageToken) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            return yield getProjectPage("/v1beta1/availableProjects", {
                responseKey: "projectInfo",
                pageSize,
                pageToken,
            });
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError("Failed to list available Google Cloud Platform projects. See firebase-debug.log for more info.", { exit: 2, original: err });
        }
    });
}
exports.getAvailableCloudProjectPage = getAvailableCloudProjectPage;
function listFirebaseProjects(pageSize) {
    return __awaiter(this, void 0, void 0, function* () {
        const projects = [];
        let nextPageToken;
        do {
            const projectPage = yield getFirebaseProjectPage(pageSize, nextPageToken);
            projects.push(...projectPage.projects);
            nextPageToken = projectPage.nextPageToken;
        } while (nextPageToken);
        return projects;
    });
}
exports.listFirebaseProjects = listFirebaseProjects;
function getFirebaseProject(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const response = yield api.request("GET", `/v1beta1/projects/${projectId}`, {
                auth: true,
                origin: api.firebaseApiOrigin,
                timeout: TIMEOUT_MILLIS,
            });
            return response.body;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError(`Failed to get Firebase project ${projectId}. ` +
                "Please make sure the project exists and your account has permission to access it.", { exit: 2, original: err });
        }
    });
}
exports.getFirebaseProject = getFirebaseProject;
