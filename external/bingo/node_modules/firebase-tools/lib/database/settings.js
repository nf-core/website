"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const defaultWriteSizeLimit = {
    path: "defaultWriteSizeLimit",
    description: `
      Set a limit for the size of each write request: small, medium, large or unlimited.
      If you choose 'unlimited', any and all write requests will be allowed, potentially
      blocking subsequent write requests while the database processes any large write
      requests. For example, deleting data at the database's root
      can't be reverted and the database will be unavailable until the delete is finished.
      To avoid blocking large write requests without running the risk of hanging your
      database, you can set this limit to small (target=10s), medium (target=30s), large (target=60s).
      Realtime Database estimates the size of each write request and aborts
      requests that will take longer than the target time.
  `,
    parseInput: (input) => {
        switch (input) {
            case "small":
            case "medium":
            case "large":
            case "unlimited":
                return input;
            default:
                return undefined;
        }
    },
    parseInputErrorMessge: "defaultWriteSizeLimit must be either small, medium, large or unlimited. (tiny is not allowed)",
};
const strictTriggerValidation = {
    path: "strictTriggerValidation",
    description: `
      Strict validation is enabled by default for write operations that trigger
      events. Any write operations that trigger more than 1000 Cloud Functions or a
      single event greater than 1 MB in size will fail and return an error reporting
      the limit that was hit. This might mean that some Cloud Functions aren't
      triggered at all if they fail the pre-validation.

      If you're performing a larger write operation (for example, deleting your
      entire database), you might want to disable this validation, as the errors
      themselves might block the operation.
  `,
    parseInput: (input) => {
        switch (input) {
            case "true":
                return true;
            case "false":
                return false;
            default:
                return undefined;
        }
    },
    parseInputErrorMessge: "strictTriggerValidation must be 'true' or 'false'",
};
exports.DATABASE_SETTINGS = new Map();
exports.DATABASE_SETTINGS.set(defaultWriteSizeLimit.path, defaultWriteSizeLimit);
exports.DATABASE_SETTINGS.set(strictTriggerValidation.path, strictTriggerValidation);
exports.HELP_TEXT = "\nAvailable Settings:\n" +
    Array.from(exports.DATABASE_SETTINGS.values())
        .map((setting) => `  ${setting.path}:${setting.description}`)
        .join("");
exports.INVALID_PATH_ERROR = `Path must be one of ${Array.from(exports.DATABASE_SETTINGS.keys()).join(", ")}.`;
