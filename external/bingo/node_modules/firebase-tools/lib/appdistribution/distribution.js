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
const fs = require("fs-extra");
const error_1 = require("../error");
const crypto = require("crypto");
var DistributionFileType;
(function (DistributionFileType) {
    DistributionFileType["IPA"] = "ipa";
    DistributionFileType["APK"] = "apk";
})(DistributionFileType = exports.DistributionFileType || (exports.DistributionFileType = {}));
class Distribution {
    constructor(path) {
        this.path = path;
        if (!path) {
            throw new error_1.FirebaseError("must specify a distribution file");
        }
        const distributionType = path.split(".").pop();
        if (distributionType !== DistributionFileType.IPA &&
            distributionType !== DistributionFileType.APK) {
            throw new error_1.FirebaseError("unsupported distribution file format, should be .ipa or .apk");
        }
        if (!fs.existsSync(path)) {
            throw new error_1.FirebaseError(`File ${path} does not exist: verify that file points to a distribution`);
        }
        this.path = path;
        this.fileType = distributionType;
    }
    fileSize() {
        return fs.statSync(this.path).size;
    }
    readStream() {
        return fs.createReadStream(this.path);
    }
    platform() {
        switch (this.fileType) {
            case DistributionFileType.IPA:
                return "ios";
            case DistributionFileType.APK:
                return "android";
            default:
                throw new error_1.FirebaseError("Unsupported distribution file format, should be .ipa or .apk");
        }
    }
    releaseHash() {
        return __awaiter(this, void 0, void 0, function* () {
            return new Promise((resolve) => {
                const hash = crypto.createHash("sha1");
                const stream = this.readStream();
                stream.on("data", (data) => hash.update(data));
                stream.on("end", () => {
                    return resolve(hash.digest("hex"));
                });
            });
        });
    }
}
exports.Distribution = Distribution;
