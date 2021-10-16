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
const inquirer = require("inquirer");
const _ = require("lodash");
const error_1 = require("./error");
function prompt(options, questions) {
    return __awaiter(this, void 0, void 0, function* () {
        const prompts = [];
        for (const question of questions) {
            if (question.name && options[question.name] === undefined) {
                prompts.push(question);
            }
        }
        if (prompts.length && options.nonInteractive) {
            const missingOptions = _.uniq(_.map(prompts, "name")).join(", ");
            throw new error_1.FirebaseError(`Missing required options (${missingOptions}) while running in non-interactive mode`, {
                children: prompts,
                exit: 1,
            });
        }
        const answers = yield inquirer.prompt(prompts);
        _.forEach(answers, (v, k) => {
            options[k] = v;
        });
        return answers;
    });
}
exports.prompt = prompt;
function promptOnce(question) {
    return __awaiter(this, void 0, void 0, function* () {
        question.name = question.name || "question";
        const answers = yield prompt({}, [question]);
        return answers[question.name];
    });
}
exports.promptOnce = promptOnce;
