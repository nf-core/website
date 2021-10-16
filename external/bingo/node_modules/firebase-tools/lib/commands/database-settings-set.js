"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const request = require("request");
const responseToError = require("../responseToError");
const command_1 = require("../command");
const error_1 = require("../error");
const requirePermissions_1 = require("../requirePermissions");
const utils = require("../utils");
const api = require("../api");
const requireInstance = require("../requireInstance");
const settings_1 = require("../database/settings");
exports.default = new command_1.Command("database:settings:set <path> <value>")
    .description("set the realtime database setting at path.")
    .option("--instance <instance>", "use the database <instance>.firebaseio.com (if omitted, use default database instance)")
    .help(settings_1.HELP_TEXT)
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.update"])
    .before(requireInstance)
    .action((path, value, options) => {
    const setting = settings_1.DATABASE_SETTINGS.get(path);
    if (setting === undefined) {
        return utils.reject(settings_1.INVALID_PATH_ERROR, { exit: 1 });
    }
    const parsedValue = setting.parseInput(value);
    if (parsedValue === undefined) {
        return utils.reject(setting.parseInputErrorMessge, { exit: 1 });
    }
    return new Promise((resolve, reject) => {
        const url = utils.addSubdomain(api.realtimeOrigin, options.instance) + "/.settings/" + path + ".json";
        const reqOptions = {
            url,
            body: parsedValue,
            json: true,
        };
        return api.addRequestHeaders(reqOptions).then((reqOptionsWithToken) => {
            request.put(reqOptionsWithToken, (err, res, body) => {
                if (err) {
                    return reject(new error_1.FirebaseError(`Unexpected error fetching configs at ${path}`, {
                        exit: 2,
                        original: err,
                    }));
                }
                else if (res.statusCode >= 400) {
                    return reject(responseToError(res, body));
                }
                utils.logSuccess("Successfully set setting.");
                utils.logSuccess(`For database instance ${options.instance}\n\t ${path} = ${JSON.stringify(parsedValue)}`);
                resolve();
            });
        });
    });
});
