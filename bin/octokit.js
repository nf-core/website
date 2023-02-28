import * as dotenv from 'dotenv';
import { Octokit } from 'octokit';

// see https://github.com/motdotla/dotenv#how-do-i-use-dotenv-with-import
dotenv.config();

const octokit = new Octokit({
  auth: process.env.GITHUB_TOKEN,
});
export default octokit;
