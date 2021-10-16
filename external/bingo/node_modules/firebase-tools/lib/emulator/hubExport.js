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
const path = require("path");
const fs = require("fs");
const api = require("../api");
const types_1 = require("./types");
const registry_1 = require("./registry");
const error_1 = require("../error");
const hub_1 = require("./hub");
const downloadableEmulators_1 = require("./downloadableEmulators");
class HubExport {
    constructor(projectId, exportPath) {
        this.projectId = projectId;
        this.exportPath = exportPath;
    }
    static readMetadata(exportPath) {
        const metadataPath = path.join(exportPath, this.METADATA_FILE_NAME);
        if (!fs.existsSync(metadataPath)) {
            return undefined;
        }
        return JSON.parse(fs.readFileSync(metadataPath, "utf8").toString());
    }
    exportAll() {
        return __awaiter(this, void 0, void 0, function* () {
            const toExport = types_1.ALL_EMULATORS.filter(this.shouldExport);
            if (toExport.length === 0) {
                throw new error_1.FirebaseError("No running emulators support import/export.");
            }
            const metadata = {
                version: hub_1.EmulatorHub.CLI_VERSION,
            };
            if (this.shouldExport(types_1.Emulators.FIRESTORE)) {
                metadata.firestore = {
                    version: downloadableEmulators_1.getDownloadDetails(types_1.Emulators.FIRESTORE).version,
                    path: "firestore_export",
                    metadata_file: "firestore_export/firestore_export.overall_export_metadata",
                };
                yield this.exportFirestore(metadata);
            }
            const metadataPath = path.join(this.exportPath, HubExport.METADATA_FILE_NAME);
            fs.writeFileSync(metadataPath, JSON.stringify(metadata));
        });
    }
    exportFirestore(metadata) {
        return __awaiter(this, void 0, void 0, function* () {
            const firestoreInfo = registry_1.EmulatorRegistry.get(types_1.Emulators.FIRESTORE).getInfo();
            const firestoreHost = `http://${firestoreInfo.host}:${firestoreInfo.port}`;
            const firestoreExportBody = {
                database: `projects/${this.projectId}/databases/(default)`,
                export_directory: this.exportPath,
                export_name: metadata.firestore.path,
            };
            return api.request("POST", `/emulator/v1/projects/${this.projectId}:export`, {
                origin: firestoreHost,
                json: true,
                data: firestoreExportBody,
            });
        });
    }
    shouldExport(e) {
        return types_1.IMPORT_EXPORT_EMULATORS.indexOf(e) >= 0 && registry_1.EmulatorRegistry.isRunning(e);
    }
}
exports.HubExport = HubExport;
HubExport.METADATA_FILE_NAME = "firebase-export-metadata.json";
