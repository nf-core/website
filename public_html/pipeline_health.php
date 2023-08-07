<?php
$title = 'Repository health';
$subtitle = 'Check GitHub settings for all sanger-tol Nextflow pipelines repositories';
$mainpage_container = false;
include '../includes/header.php';

// Refresh cache?
function is_refresh_cache($repo = null, $any_repo = false) {
    if (!isset($_GET['action']) || $_GET['action'] != 'refresh') {
        return false;
    }
    if ($any_repo || (isset($_GET['repos']) && $_GET['repos'] == 'all')) {
        return true;
    }
    if ($repo && isset($_GET['repos']) && $_GET['repos'] == $repo) {
        return true;
    }
    return false;
}

// Fix repo?
function is_fix_repo($repo = null) {
    if (!isset($_GET['action']) || $_GET['action'] != 'fix') {
        return false;
    }
    if (isset($_GET['repos']) && $_GET['repos'] == 'all') {
        return true;
    }
    if ($repo && isset($_GET['repos']) && $_GET['repos'] == $repo) {
        return true;
    }
    return false;
}

// Get auth secrets
$config = parse_ini_file('../config.ini');
define('GH_AUTH', base64_encode($config['github_username'] . ':' . $config['github_access_token']));

// Load pipelines JSON
$pipelines_json = json_decode(file_get_contents('pipelines.json'))->remote_workflows;

// Placeholders
$pipelines = [];
$core_repos = [];

// HTTP header to use on GitHub API GET requests
define(
    'GH_API_OPTS',
    stream_context_create([
        'http' => [
            'method' => 'GET',
            'header' => [
                'User-Agent: PHP',
                'Accept:application/vnd.github.mercy-preview+json', // Needed to get topics (keywords) for now
                'Accept:application/vnd.github.luke-cage-preview+json', // Needed to get protected branch required reviews
                'Authorization: Basic ' . GH_AUTH,
            ],
        ],
    ]),
);

// Base repo health class
class RepoHealth {
    // Init vars
    public $name;
    public $refresh = false;
    public $cache_base;
    public function __construct($name) {
        $this->name = $name;
        $this->refresh = is_refresh_cache($this->name);

        // Cache filenames
        $this->cache_base = dirname(dirname(__FILE__)) . '/api_cache/pipeline_health';
        $this->gh_repo_cache = $this->cache_base . '/repo_' . $this->name . '.json';
        $this->gh_release_cache = $this->cache_base . '/release_' . $this->name . '.json';
        $this->gh_all_branches_cache = $this->cache_base . '/branches_' . $this->name . '.json';
    }

    // Names of required CI checks. These are added to whatever already exists.
    public $required_status_check_contexts = [
        'Prettier',
        'EditorConfig',
        'nf-core',
        'Run pipeline with test data',
        // NOTE - doesn't seem to be any way to get the "available" contexts through GitHub API
        // If we really want to do this, might have to query the repo contents..??
    ];

    // Names of old CI tests that must not be present any more
    public $required_remove_status_check_contexts = [
        'continuous-integration/travis-ci',
        'test',
        'YAML',
        'Markdown',
        'Run workflow tests',
    ];
    public $branch_exist_tests = ['main'];
    public $branches_protection = ['main'];
    // public $branch_template_protection = false;
    public $branch_default = 'main';
    public $required_topics = ['nextflow'];
    public $web_url = 'https://pipelines.tol.sanger.ac.uk';
    public $test_names;
    public $test_descriptions;
    public $test_urls;

    // Data vars
    public $gh_repo;
    public $gh_release;
    public $gh_teams = [];
    public $gh_branches;
    public $gh_branch_main;
    public $gh_branch_dev;
    public $gh_branch_TEMPLATE;

    // Test result variables
    public $repo_wikis;
    public $repo_issues;
    public $repo_merge_commits;
    public $repo_merge_rebase;
    public $repo_merge_squash;
    public $repo_default_branch;
    public $repo_keywords;
    public $repo_description;
    public $repo_url;
    public $team_nextflow_all;
    public $team_nextflow_admin;

    // Branch test vars
    public $branch_main_exists;
    public $branch_dev_exists;
    public $branch_template_exists;
    public $branch_main_strict_updates;
    public $branch_main_required_ci;
    public $branch_main_stale_reviews;
    public $branch_main_code_owner_reviews;
    public $branch_main_required_num_reviews;
    public $branch_main_enforce_admins;
    public $branch_dev_strict_updates;
    public $branch_dev_required_ci;
    public $branch_dev_stale_reviews;
    public $branch_dev_code_owner_reviews;
    public $branch_dev_required_num_reviews;
    public $branch_dev_enforce_admins;
    // public $branch_template_restrict_push;

    public function get_data() {
        $this->get_repo_data();
        $this->get_branch_data();
    }
    public function run_tests() {
        $this->test_repo();
        $this->test_teams();
        $this->test_branch_exists();
    }
    public function fix_tests() {
        if (is_fix_repo($this->name)) {
            $this->fix_repo();
            $this->fix_topics();
            $this->fix_teams();
            // Done! Refresh the test statuses
            $this->run_tests();
        }
    }

    public function get_repo_data() {
        // Super annoyingly, the teams call misses just one or two keys we need :(
        if (is_null($this->gh_repo) || !isset($this->gh_repo->allow_merge_commit)) {
            if (file_exists($this->gh_repo_cache) && !$this->refresh) {
                $this->gh_repo = json_decode(file_get_contents($this->gh_repo_cache));
            } else {
                $gh_repo_url = 'https://api.github.com/repos/sanger-tol/' . $this->name;
                $this->gh_repo = json_decode(file_get_contents($gh_repo_url, false, GH_API_OPTS));
                $this->_save_cache_data($this->gh_repo_cache, $this->gh_repo);
            }
        }
    }

    public function get_release_data() {
        // Currently only used to get last release for tools, as have others from pipelines.json
        if (file_exists($this->gh_release_cache) && !$this->refresh) {
            $this->gh_release = json_decode(file_get_contents($this->gh_release_cache));
        } else {
            $gh_release_url = 'https://api.github.com/repos/nf-core/' . $this->name . '/releases/latest';
            $this->gh_release = json_decode(file_get_contents($gh_release_url, false, GH_API_OPTS));
            $this->_save_cache_data($this->gh_release_cache, $this->gh_release);
        }
    }

    public function get_branch_data() {
        // List of all branches
        if (file_exists($this->gh_all_branches_cache) && !$this->refresh) {
            $this->gh_branches = json_decode(file_get_contents($this->gh_all_branches_cache));
        } else {
            $gh_branch_url = 'https://api.github.com/repos/sanger-tol/' . $this->name . '/branches';
            $this->gh_branches = json_decode(@file_get_contents($gh_branch_url, false, GH_API_OPTS));
            $this->_save_cache_data($this->gh_all_branches_cache, $this->gh_branches);
        }

        // Details of branch protection for main and dev
        foreach (['main', 'dev', 'TEMPLATE'] as $branch) {
            $gh_branch_cache = $this->cache_base . '/branch_' . $this->name . '_' . $branch . '.json';
            if (file_exists($gh_branch_cache) && !$this->refresh) {
                $gh_branch = json_decode(file_get_contents($gh_branch_cache));
                if (is_object($gh_branch)) {
                    $this->{'gh_branch_' . $branch} = $gh_branch;
                } else {
                }
            } else {
                $gh_branch_url =
                    'https://api.github.com/repos/sanger-tol/' . $this->name . '/branches/' . $branch . '/protection';
                $gh_branch = json_decode(@file_get_contents($gh_branch_url, false, GH_API_OPTS));
                if (strpos($http_response_header[0], 'HTTP/1.1 200') !== false && is_object($gh_branch)) {
                    $this->{'gh_branch_' . $branch} = $gh_branch;
                    $this->_save_cache_data($gh_branch_cache, $this->{'gh_branch_' . $branch});
                } else {
                    // Write an empty cache file
                    if (strpos($http_response_header[0], 'HTTP/1.1 404') === false) {
                        // A 404 is fine, that just means that there is no branch protection. Warn if anything else.
                        echo '<div class="alert alert-danger">Could not fetch branch protection data for <code>' .
                            $this->name .
                            '</code> - <code>' .
                            $branch .
                            '</code><pre>' .
                            print_r($http_response_header, true) .
                            '</pre><pre>' .
                            print_r($gh_branch, true) .
                            '</pre></div>';
                    }
                    $this->_save_cache_data($gh_branch_cache, '');
                }
            }
        }
    }

    public function test_topics() {
        $topics_pass = true;
        foreach ($this->required_topics as $top) {
            if (!in_array($top, $this->gh_repo->topics)) {
                $topics_pass = false;
                break;
            }
        }
        return $topics_pass;
    }
    public function test_repo() {
        if (isset($this->gh_repo->has_wiki)) {
            $this->repo_wikis = !$this->gh_repo->has_wiki;
        }
        if (isset($this->gh_repo->has_issues)) {
            $this->repo_issues = $this->gh_repo->has_issues;
        }
        if (isset($this->gh_repo->allow_merge_commit)) {
            $this->repo_merge_commits = $this->gh_repo->allow_merge_commit;
        }
        if (isset($this->gh_repo->allow_rebase_merge)) {
            $this->repo_merge_rebase = $this->gh_repo->allow_rebase_merge;
        }
        if (isset($this->gh_repo->allow_squash_merge)) {
            $this->repo_merge_squash = $this->gh_repo->allow_squash_merge;
        }
        if (isset($this->gh_repo->default_branch)) {
            $this->repo_default_branch = $this->gh_repo->default_branch == $this->branch_default;
        }
        if (isset($this->gh_repo->topics)) {
            $this->repo_keywords = $this->test_topics();
        }
        if (isset($this->gh_repo->description)) {
            $this->repo_description = $this->gh_repo->description;
        }
        if (isset($this->gh_repo->homepage)) {
            $this->repo_url = $this->gh_repo->homepage == $this->web_url;
        }
    }
    public function test_teams() {
        $this->team_nextflow_all = isset($this->gh_teams['nextflow_all']) ? $this->gh_teams['nextflow_all']->push : false;
        $this->team_nextflow_admin = isset($this->gh_teams['nextflow_admin']) ? $this->gh_teams['nextflow_admin']->admin : false;
    }
    public function test_branch_exists() {
        // Check that branches exist
        if (isset($this->gh_branches)) {
            $this->branch_main_exists = false;
            $this->branch_dev_exists = false;
            $this->branch_template_exists = false;
            foreach ($this->gh_branches as $branch) {
                if (in_array(strtolower($branch->name), $this->branch_exist_tests)) {
                    $this->{'branch_' . strtolower($branch->name) . '_exists'} = true;
                }
            }
            // Check that the branch that should be default actually exists
            if (!$this->{'branch_' . strtolower($this->branch_default) . '_exists'}) {
                $this->repo_default_branch = -1;
            }
        }
    }
    public function test_branch_protection() {
        // Test branch protection for main and dev
        foreach ($this->branches_protection as $branch) {
            $prs_required = $branch == 'main' ? 2 : 1;
            if (!$this->{'branch_' . $branch . '_exists'}) {
                $this->{'branch_' . $branch . '_strict_updates'} = -1;
                $this->{'branch_' . $branch . '_required_ci'} = -1;
                $this->{'branch_' . $branch . '_stale_reviews'} = -1;
                $this->{'branch_' . $branch . '_code_owner_reviews'} = -1;
                $this->{'branch_' . $branch . '_required_num_reviews'} = -1;
                $this->{'branch_' . $branch . '_enforce_admins'} = -1;
                $this->test_descriptions['branch_' . $branch . '_strict_updates'] = $branch . ' branch does not exist';
                $this->test_descriptions['branch_' . $branch . '_required_ci'] = $branch . ' branch does not exist';
                $this->test_descriptions['branch_' . $branch . '_stale_reviews'] = $branch . ' branch does not exist';
                $this->test_descriptions['branch_' . $branch . '_code_owner_reviews'] =
                    $branch . ' branch does not exist';
                $this->test_descriptions['branch_' . $branch . '_required_num_reviews'] =
                    $branch . ' branch does not exist';
                $this->test_descriptions['branch_' . $branch . '_enforce_admins'] = $branch . ' branch does not exist';
                continue;
            }
            $data = $this->{'gh_branch_' . $branch};

            if (!isset($data->required_status_checks)) {
                $this->{'branch_' . $branch . '_strict_updates'} = false;
                $this->{'branch_' . $branch . '_required_ci'} = false;
            } else {
                // Don't require branches to be up to date before merging.
                $this->{'branch_' . $branch . '_strict_updates'} = $data->required_status_checks->strict == false;
                // At least the minimum set of required CI tests
                $has_min_ci_checks = !array_diff(
                    $this->required_status_check_contexts,
                    $data->required_status_checks->contexts,
                );
                // None of the required tests that we don't want
                $has_removed_ci_checks =
                    count(
                        array_diff(
                            $this->required_remove_status_check_contexts,
                            $data->required_status_checks->contexts,
                        ),
                    ) == count($this->required_remove_status_check_contexts);
                // Only pass if both of the above pass
                $this->{'branch_' . $branch . '_required_ci'} = $has_min_ci_checks && $has_removed_ci_checks;
            }
            if (!isset($data->required_pull_request_reviews)) {
                $this->{'branch_' . $branch . '_stale_reviews'} = false;
                $this->{'branch_' . $branch . '_code_owner_reviews'} = false;
                $this->{'branch_' . $branch . '_required_num_reviews'} = false;
            } else {
                // Don't mark reviews as stale on new commits
                $this->{'branch_' . $branch . '_stale_reviews'} =
                    $data->required_pull_request_reviews->dismiss_stale_reviews == false;
                // Don't require reviews from code owners
                $this->{'branch_' . $branch . '_code_owner_reviews'} =
                    $data->required_pull_request_reviews->require_code_owner_reviews == false;
                // Require 1 or 2 reviews
                $this->{'branch_' . $branch . '_required_num_reviews'} =
                    $data->required_pull_request_reviews->required_approving_review_count == $prs_required;
            }
            // Don't include administrators
            if (!isset($data->enforce_admins)) {
                $this->{'branch_' . $branch . '_enforce_admins'} = false;
            } else {
                $this->{'branch_' . $branch . '_enforce_admins'} = $data->enforce_admins->enabled == false;
            }
        }
        // Tests specifically for the TEMPLATE branch
        /*
        if (!$this->branch_template_exists) {
            $this->branch_template_restrict_push = -1;
            $this->test_descriptions['branch_template_restrict_push'] = 'TEMPLATE branch does not exist';
        } else {
            $data = $this->gh_branch_TEMPLATE;
            $this->branch_template_restrict_push = false;
            if (isset($data->restrictions)) {
                if (count($data->restrictions->users) == 1) {
                    if ($data->restrictions->users[0]->login == 'nf-core-bot') {
                        $this->branch_template_restrict_push = true;
                    }
                }
            }
        }
        */
    }

    private function fix_repo() {
        // https://developer.github.com/v3/repos/#edit
        $payload = [];
        if (!$this->repo_wikis) {
            $payload['has_wiki'] = false;
        }
        if (!$this->repo_issues) {
            $payload['has_issues'] = true;
        }
        if (!$this->repo_merge_commits) {
            $payload['allow_merge_commit'] = true;
        }
        if (!$this->repo_merge_rebase) {
            $payload['allow_rebase_merge'] = true;
        }
        if (!$this->repo_merge_squash) {
            $payload['allow_squash_merge'] = true;
        }
        if (!$this->repo_default_branch) {
            $payload['default_branch'] = $this->branch_default;
        }
        if (!$this->repo_url) {
            $payload['homepage'] = $this->web_url;
        }
        if (count($payload) > 0) {
            $gh_edit_repo_url = 'https://api.github.com/repos/sanger-tol/' . $this->name;
            $updated_data = $this->_send_gh_api_data($gh_edit_repo_url, $payload, 'PATCH');
            if ($updated_data) {
                $this->gh_repo = $updated_data;
                $this->_save_cache_data($this->gh_repo_cache, $this->gh_repo);
            }
        }
    }

    private function fix_topics() {
        // https://developer.github.com/v3/repos/#replace-all-topics-for-a-repository
        if (!$this->repo_keywords) {
            $topics = [
                'names' => array_values(array_unique(array_merge($this->gh_repo->topics, $this->required_topics))),
            ];
            $gh_edit_topics_url = 'https://api.github.com/repos/sanger-tol/' . $this->name . '/topics';
            $updated_data = $this->_send_gh_api_data($gh_edit_topics_url, $topics, 'PUT');
            if ($updated_data) {
                $this->gh_repo->topics = $updated_data->names;
                $this->_save_cache_data($this->gh_repo_cache, $this->gh_repo);
            }
        }
    }

    private function fix_teams() {
        $this->fix_team('nextflow_all');
        $this->fix_team('nextflow_admin');
    }
    private function fix_team($team) {
        global $gh_team_ids;
        global $updated_teams;
        if (!$this->{'team_' . $team}) {
            $payload = [];
            if ($team == 'nextflow_admin') {
                $payload = ['permission' => 'admin'];
            }
            if ($team == 'nextflow_all') {
                $payload = ['permission' => 'push'];
            }
            $gh_edit_team_url = 'https://api.github.com/teams/' . $gh_team_ids[$team] . '/repos/sanger-tol/' . $this->name;
            if ($this->_send_gh_api_data($gh_edit_team_url, $payload, 'PUT')) {
                $updated_teams[$team] = true;
            }
        }
    }

    public function fix_branch_protection() {
        // Fix branch protection for main and dev
        foreach ($this->branches_protection as $branch) {
            // Convenience vars for test results
            $test_results = [
                $this->{'branch_' . $branch . '_enforce_admins'},
                $this->{'branch_' . $branch . '_strict_updates'},
                $this->{'branch_' . $branch . '_required_ci'},
                $this->{'branch_' . $branch . '_stale_reviews'},
                $this->{'branch_' . $branch . '_code_owner_reviews'},
                $this->{'branch_' . $branch . '_required_num_reviews'},
            ];

            // Only run if we have at least one test failure
            if (count(array_keys($test_results, true)) != count($test_results)) {
                // Add needed required-CI tests to what's already there if we have something
                if (
                    is_object($this->{'gh_branch_' . $branch}) &&
                    isset($this->{'gh_branch_' . $branch}->required_status_checks)
                ) {
                    $contexts = array_values(
                        array_unique(
                            array_merge(
                                $this->required_status_check_contexts,
                                $this->{'gh_branch_' . $branch}->required_status_checks->contexts,
                            ),
                        ),
                    );
                } else {
                    $contexts = $this->required_status_check_contexts;
                }
                // Remove any old contexts that we don't want
                $contexts = array_diff($contexts, $this->required_remove_status_check_contexts);

                $payload = [
                    'enforce_admins' => false,
                    'required_status_checks' => [
                        'strict' => false,
                        'contexts' => array_values($contexts),
                    ],
                    'required_pull_request_reviews' => [
                        'dismiss_stale_reviews' => false,
                        'require_code_owner_reviews' => false,
                        'required_approving_review_count' => $branch == 'main' ? 2 : 1,
                    ],
                    'restrictions' => null,
                ];
                $gh_edit_branch_protection_url =
                    'https://api.github.com/repos/sanger-tol/' . $this->name . '/branches/' . $branch . '/protection';
                $updated_data = $this->_send_gh_api_data($gh_edit_branch_protection_url, $payload, 'PUT');
                if ($updated_data) {
                    $this->{'gh_branch_' . $branch} = $updated_data;
                    $gh_branch_cache = $this->cache_base . '/branch_' . $this->name . '_' . $branch . '.json';
                    $this->_save_cache_data($gh_branch_cache, $this->{'gh_branch_' . $branch});
                }
            }
        }

        // Fix TEMPLATE branch protection
        /*
        if ($this->branch_template_protection) {
            // Only run if the test failed
            if ($this->branch_template_restrict_push === false) {
                $payload = [
                    'enforce_admins' => false,
                    'required_status_checks' => null,
                    'required_pull_request_reviews' => null,
                    'restrictions' => [
                        'users' => ['nf-core-bot'],
                        'teams' => [],
                    ],
                ];
                // Push to GitHub API
                $gh_template_branch_protection_url =
                    'https://api.github.com/repos/sanger-tol/' . $this->name . '/branches/TEMPLATE/protection';
                $updated_data = $this->_send_gh_api_data($gh_template_branch_protection_url, $payload, 'PUT');
                if ($updated_data) {
                    $this->gh_branch_TEMPLATE = $updated_data;
                    $gh_branch_cache = $this->cache_base . '/branch_' . $this->name . '_TEMPLATE.json';
                    $this->_save_cache_data($gh_branch_cache, $this->gh_branch_TEMPLATE);
                }
            }
        }
        */
    }

    public function _send_gh_api_data($url, $content, $method = 'POST') {
        $context = stream_context_create([
            'http' => [
                'method' => $method,
                'header' => [
                    'Content-Type: application/json',
                    'User-Agent: PHP',
                    'Accept:application/vnd.github.mercy-preview+json', // Needed to get topics (keywords) for now
                    'Accept:application/vnd.github.luke-cage-preview+json', // Needed to get protected branch required reviews
                    'Authorization: Basic ' . GH_AUTH,
                ],
                'content' => json_encode($content),
            ],
        ]);
        $result = json_decode(file_get_contents($url, false, $context));
        if (strpos($http_response_header[0], 'HTTP/1.1 204') !== false) {
            return true;
        } elseif (strpos($http_response_header[0], 'HTTP/1.1 200') !== false) {
            return $result;
        } else {
            echo '<div class="alert alert-danger m-3">
        <strong class="me-2">Error with GitHub API</strong>
        There was a problem with the following URL: <code class="me-2">' .
                $url .
                '</code> (<code>' .
                $method .
                '</code>)
        <details>
          <pre>' .
                json_encode($http_response_header, JSON_PRETTY_PRINT) .
                '</pre>
          <pre>' .
                json_encode($content, JSON_PRETTY_PRINT) .
                '</pre>
        </details>
      </div>';
            return false;
        }
    }

    public function _save_cache_data($path, $data, $encode_json = true) {
        if (!file_exists(dirname($path))) {
            mkdir(dirname($path), 0777, true);
        }
        if ($encode_json) {
            $data_json = json_encode($data, JSON_PRETTY_PRINT) . "\n";
        } else {
            $data_json = $data;
        }
        file_put_contents($path, $data_json);
    }

    public function print_table_cell($test_name) {
        $test_url = $this->test_urls[$test_name];
        $test_url = str_replace('{repo}', $this->name, $test_url);
        if (property_exists($this, 'last_release') && $this->last_release) {
            $test_url = str_replace('{latest-tag}', $this->last_release->tag_name, $test_url);
        }
        if (is_null($this->$test_name)) {
            echo '<td class="table-secondary text-center" title="' .
                $this->name .
                ': ' .
                $this->test_descriptions[$test_name] .
                '" data-bs-toggle="tooltip" data-html="true">
        <a href="' .
                $test_url .
                '" class="d-block" target="_blank"><i class="fas fa-question text-secondary"></i></a>
      </td>';
        } elseif ($this->$test_name === -1) {
            echo '<td class="table-success text-center" title="' .
                $this->name .
                ': ' .
                $this->test_descriptions[$test_name] .
                '" data-bs-toggle="tooltip" data-html="true">
        <a href="' .
                $test_url .
                '" class="d-block text-secondary text-decoration-none" target="_blank">&mdash;</a>
      </td>';
        } elseif ($this->$test_name) {
            echo '<td class="table-success text-center" title="' .
                $this->name .
                ': ' .
                $this->test_descriptions[$test_name] .
                '" data-bs-toggle="tooltip" data-html="true">
        <a href="' .
                $test_url .
                '" class="d-block" target="_blank"><i class="fas fa-check text-success"></i></a>
      </td>';
        } else {
            echo '<td class="table-danger text-center" title="' .
                $this->name .
                ': ' .
                $this->test_descriptions[$test_name] .
                '" data-bs-toggle="tooltip" data-html="true">
        <a href="' .
                $test_url .
                '" class="d-block" target="_blank"><i class="fas fa-times text-danger"></i></a>
      </td>';
        }
    }
}

// Pipeline health class
class PipelineHealth extends RepoHealth {
    // URL should point to pipeline page
    public function __construct($name) {
        parent::__construct($name);
        $this->web_url = 'https://pipelines.tol.sanger.ac.uk/' . $this->name;
    }
    // We need more branches in pipelines
    public $branch_exist_tests = ['template', 'dev', 'main']; // lower case
    public $branches_protection = ['dev', 'main'];
    // public $branch_template_protection = true;
    // Keywords should also include nextflow, workflow and pipeline
    public $required_topics = ['nextflow', 'pipeline'];
    // JSON Schema / DSL2 modules directory
    public $has_json_schema;
    public $has_dsl2_modules_dir;
    // Variables for release tests
    public $has_release;
    public $last_release;
    public $release_after_tools;
    public $main_is_release;

    // Extra pipeline-specific tests
    public function run_tests() {
        parent::run_tests();
        $this->test_branch_protection();
        $this->test_files_exist();
        $this->test_releases();
    }

    public function check_url($url) {
        $headers = get_headers($url);
        return substr($headers[0], 9, 3) == '200';
    }
    public function test_files_exist() {
        // No releases - always check the dev branch (no caching)
        if (!$this->has_release) {
            $this->has_json_schema = $this->check_url(
                'https://raw.githubusercontent.com/sanger-tol/' . $this->name . '/dev/nextflow_schema.json',
            );
            $this->has_dsl2_modules_dir = $this->check_url(
                'https://github.com/sanger-tol/' . $this->name . '/tree/dev/modules',
            );
        }

        // Check last release, with caching
        else {
            $check_404_cache =
                $this->cache_base . '/files_404_' . $this->name . '_' . $this->last_release->tag_name . '.json';
            // Load cache
            if (file_exists($check_404_cache) && !$this->refresh) {
                $files_404_cache = json_decode(file_get_contents($check_404_cache));
                $this->has_json_schema = $files_404_cache->json_schema;
                $this->has_dsl2_modules_dir = $files_404_cache->dsl2_modules_dir;
            } else {
                // Check if the files exist
                $this->has_json_schema = $this->check_url(
                    'https://raw.githubusercontent.com/sanger-tol/' .
                        $this->name .
                        '/' .
                        $this->last_release->tag_name .
                        '/nextflow_schema.json',
                );
                $this->has_dsl2_modules_dir = $this->check_url(
                    'https://github.com/sanger-tol/' . $this->name . '/tree/' . $this->last_release->tag_name . '/modules',
                );
                // Save the cache
                $files_404_cache = [
                    'json_schema' => $this->has_json_schema,
                    'dsl2_modules_dir' => $this->has_dsl2_modules_dir,
                ];
                $this->_save_cache_data($check_404_cache, $files_404_cache);
            }
        }

        if (file_exists($this->gh_all_branches_cache) && !$this->refresh) {
            $this->gh_branches = json_decode(file_get_contents($this->gh_all_branches_cache));
        } else {
            $gh_branch_url = 'https://api.github.com/repos/sanger-tol/' . $this->name . '/branches';
            $this->gh_branches = json_decode(@file_get_contents($gh_branch_url, false, GH_API_OPTS));
            $this->_save_cache_data($this->gh_all_branches_cache, $this->gh_branches);
        }
    }

    public function test_releases() {
        global $tools_last_release;
        // No releases - set to -1 and return
        if (!$this->has_release) {
            $this->release_after_tools = -1;
            $this->main_is_release = -1;
            return;
        }
        // Check if release is after last tools release
        if ($this->last_release && $tools_last_release) {
            $this->release_after_tools = strtotime($this->last_release->published_at) > strtotime($tools_last_release);
        }
        // Check if main commit hash is same as release hash
        if ($this->last_release) {
            foreach ($this->gh_branches as $branch) {
                if ($branch->name == 'main') {
                    $this->main_is_release = $this->last_release->tag_sha == $branch->commit->sha;
                }
            }
        }
    }

    // Extra pipeline-specific fixes
    public function fix_tests() {
        parent::fix_tests();
        if (is_fix_repo($this->name)) {
            $this->fix_branch_protection();
            // Done! Refresh the test statuses
            $this->run_tests();
        }
    }
}

// Core repo health class
class CoreRepoHealth extends RepoHealth {
}

// Get nf-core GitHub teams info & repos
function get_gh_team_repos($team) {
    // Globals
    global $pipelines_json;
    global $pipelines;
    global $core_repos;
    global $gh_team_ids;

    // Get team ID
    $gh_teams_cache = dirname(dirname(__FILE__)) . '/api_cache/pipeline_health/team_' . $team . '.json';
    if (file_exists($gh_teams_cache) && !is_refresh_cache(null, true)) {
        $gh_team = json_decode(file_get_contents($gh_teams_cache));
    } else {
        $gh_team_url = 'https://api.github.com/orgs/sanger-tol/teams/' . $team;
        $gh_team = json_decode(file_get_contents($gh_team_url, false, GH_API_OPTS));

        // Save for next time
        if (!file_exists(dirname($gh_teams_cache))) {
            mkdir(dirname($gh_teams_cache), 0777, true);
        }
        $gh_team_json = json_encode($gh_team, JSON_PRETTY_PRINT) . "\n";
        file_put_contents($gh_teams_cache, $gh_team_json);
    }
    $gh_team_ids[$team] = $gh_team->id;

    $gh_team_repos_cache = dirname(dirname(__FILE__)) . '/api_cache/pipeline_health/team_' . $team . '_repos.json';
    if (file_exists($gh_team_repos_cache) && !is_refresh_cache(null, true)) {
        $gh_team_repos = json_decode(file_get_contents($gh_team_repos_cache));
    } else {
        $gh_team_repos_url = 'https://api.github.com/teams/' . $gh_team->id . '/repos';
        $first_page = true;
        $next_page = false;
        $gh_team_repos = [];
        while ($first_page || $next_page) {
            // reset loop vars
            $first_page = false;
            // Get GitHub API results
            if ($next_page) {
                $gh_team_repos_url = $next_page;
            }
            $gh_team_repos = array_merge(
                $gh_team_repos,
                json_decode(file_get_contents($gh_team_repos_url, false, GH_API_OPTS)),
            );

            // Look for URL to next page of API results
            $next_page = false;
            $m_array = preg_grep('/rel="next"/', $http_response_header);
            if (count($m_array) > 0) {
                preg_match('/<([^>]+)>; rel="next"/', array_values($m_array)[0], $matches);
                if (isset($matches[1])) {
                    $next_page = $matches[1];
                }
            }
        }

        // Save for next time
        if (!file_exists(dirname($gh_team_repos_cache))) {
            mkdir(dirname($gh_team_repos_cache), 0777, true);
        }
        $gh_team_repos_json = json_encode($gh_team_repos, JSON_PRETTY_PRINT) . "\n";
        file_put_contents($gh_team_repos_cache, $gh_team_repos_json);
    }

    // Make repo health objects
    foreach ($gh_team_repos as $repo) {
        if ($repo->archived) {
            continue;
        }
        // Make a pipeline object
        $is_pipeline = false;
        foreach ($pipelines_json as $wf) {
            if ($wf->name == $repo->name) {
                if (!array_key_exists($repo->name, $pipelines)) {
                    $pipelines[$repo->name] = new PipelineHealth($repo->name);
                    $pipelines[$repo->name]->gh_repo = $repo;
                }
                $pipelines[$repo->name]->gh_teams[$team] = $repo->permissions;
                $is_pipeline = true;
            }
        }
        // Make a core repo object
        if (!$is_pipeline) {
            if (!array_key_exists($repo->name, $core_repos)) {
                $core_repos[$repo->name] = new CoreRepoHealth($repo->name);
                $core_repos[$repo->name]->gh_repo = $repo;
            }
            $core_repos[$repo->name]->gh_teams[$team] = $repo->permissions;
        }
    }
}
$gh_team_ids = [];
get_gh_team_repos('nextflow_admin');
get_gh_team_repos('nextflow_all');

// Loop through pipelines
foreach ($pipelines_json as $wf) {
    // Remove archived pipelines
    if ($wf->archived) {
        if (array_key_exists($wf->name, $pipelines)) {
            unset($pipelines[$wf->name]);
        }
    } else {
        // Add, in case there are any without team access
        if (!array_key_exists($wf->name, $pipelines)) {
            $pipelines[$wf->name] = new PipelineHealth($wf->name);
        }
        // Add data for release tests
        $pipelines[$wf->name]->has_release = false;
        if (count($wf->releases) > 0) {
            $pipelines[$wf->name]->has_release = true;
            $pipelines[$wf->name]->last_release = end($wf->releases);
        } else {
            $pipelines[$wf->name]->branch_default = 'dev';
        }
    }
}

$base_test_names = [
    'repo_wikis' => 'Wikis',
    'repo_issues' => 'Issues',
    'repo_merge_commits' => 'Merge commits',
    'repo_merge_rebase' => 'Rebase merging',
    'repo_merge_squash' => 'Squash merges',
    'repo_default_branch' => 'Default branch',
    'repo_keywords' => 'Keywords',
    'repo_description' => 'Description',
    'repo_url' => 'Repo URL',
    'team_nextflow_all' => 'Team all',
    'team_nextflow_admin' => 'Team admin',
    'branch_main_exists' => 'main: exists',
    'branch_dev_exists' => 'dev: exists',
    'branch_template_exists' => 'TEMPLATE: exists',
    'branch_main_strict_updates' => 'main: strict updates',
    'branch_main_required_ci' => 'main: required CI',
    'branch_main_stale_reviews' => 'main: stale reviews',
    'branch_main_code_owner_reviews' => 'main: code owner reviews',
    'branch_main_required_num_reviews' => 'main: 2 reviews',
    'branch_main_enforce_admins' => 'main: enforce admins',
    'branch_dev_strict_updates' => 'dev: strict updates',
    'branch_dev_required_ci' => 'dev: required CI',
    'branch_dev_stale_reviews' => 'dev: stale reviews',
    'branch_dev_code_owner_reviews' => 'dev: code owner reviews',
    'branch_dev_required_num_reviews' => 'dev: 1 review',
    'branch_dev_enforce_admins' => 'dev: enforce admins',
    // 'branch_template_restrict_push' => 'T push',
];
$base_test_descriptions = [
    'repo_wikis' => 'Disable wikis',
    'repo_issues' => 'Enable issues',
    'repo_merge_commits' => 'Allow merge commits',
    'repo_merge_rebase' => 'Allow rebase merging',
    'repo_merge_squash' => 'Allow squash merges',
    'repo_default_branch' => 'default branch main (released) or dev (no releases)',
    'repo_keywords' => 'Minimum keywords set',
    'repo_description' => 'Description must be set',
    'repo_url' => 'URL should be set to https://pipelines.tol.sanger.ac.uk',
    'team_nextflow_all' => 'Write access for sanger-tol/nextflow-all',
    'team_nextflow_admin' => 'Admin access for sanger-tol/nextflow-admin',
    'branch_main_exists' => 'main branch: branch must exist',
    'branch_dev_exists' => 'dev branch: branch must exist',
    'branch_template_exists' => 'TEMPLATE branch: branch must exist',
    'branch_main_strict_updates' => 'main branch: do not require branch to be up to date before merging',
    'branch_main_required_ci' => 'main branch: minimum set of CI tests must pass',
    'branch_main_stale_reviews' => 'main branch: reviews not marked stale after new commits',
    'branch_main_code_owner_reviews' => 'main branch: code owner reviews not required',
    'branch_main_required_num_reviews' => 'main branch: 2 reviews required',
    'branch_main_enforce_admins' => 'main branch: do not enforce rules for admins',
    'branch_dev_strict_updates' => 'dev branch: do not require branch to be up to date before merging',
    'branch_dev_required_ci' => 'dev branch: minimum set of CI tests must pass',
    'branch_dev_stale_reviews' => 'dev branch: reviews not marked stale after new commits',
    'branch_dev_code_owner_reviews' => 'dev branch: code owner reviews not required',
    'branch_dev_required_num_reviews' => 'dev branch: 1 review required',
    'branch_dev_enforce_admins' => 'dev branch: do not enforce rules for admins',
    // 'branch_template_restrict_push' => 'Restrict push to TEMPLATE to @nf-core-bot',
];
$base_test_urls = [
    'repo_wikis' => 'https://github.com/sanger-tol/{repo}/settings',
    'repo_issues' => 'https://github.com/sanger-tol/{repo}/settings',
    'repo_merge_commits' => 'https://github.com/sanger-tol/{repo}/settings',
    'repo_merge_rebase' => 'https://github.com/sanger-tol/{repo}/settings',
    'repo_merge_squash' => 'https://github.com/sanger-tol/{repo}/settings',
    'repo_default_branch' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'repo_keywords' => 'https://github.com/sanger-tol/{repo}',
    'repo_description' => 'https://github.com/sanger-tol/{repo}',
    'repo_url' => 'https://github.com/sanger-tol/{repo}',
    'team_nextflow_all' => 'https://github.com/sanger-tol/{repo}/settings/collaboration',
    'team_nextflow_admin' => 'https://github.com/sanger-tol/{repo}/settings/collaboration',
    'branch_main_exists' => 'https://github.com/sanger-tol/{repo}/branches',
    'branch_dev_exists' => 'https://github.com/sanger-tol/{repo}/branches',
    'branch_template_exists' => 'https://github.com/sanger-tol/{repo}/branches',
    'branch_main_strict_updates' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_main_required_ci' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_main_stale_reviews' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_main_code_owner_reviews' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_main_required_num_reviews' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_main_enforce_admins' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_dev_strict_updates' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_dev_required_ci' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_dev_stale_reviews' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_dev_code_owner_reviews' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_dev_required_num_reviews' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    'branch_dev_enforce_admins' => 'https://github.com/sanger-tol/{repo}/settings/branches',
    // 'branch_template_restrict_push' => 'https://github.com/sanger-tol/{repo}/settings/branches',
];
$base_merge_table_col_headings = [
    'Team access' => ['team_nextflow_all', 'team_nextflow_admin'],
    'Branches exist' => ['branch_main_exists', 'branch_dev_exists', 'branch_template_exists'],
    'Branch protection: main' => [
        'branch_main_strict_updates',
        'branch_main_required_ci',
        'branch_main_stale_reviews',
        'branch_main_code_owner_reviews',
        'branch_main_required_num_reviews',
        'branch_main_enforce_admins',
    ],
    'Branch protection: dev' => [
        'branch_dev_strict_updates',
        'branch_dev_required_ci',
        'branch_dev_stale_reviews',
        'branch_dev_code_owner_reviews',
        'branch_dev_required_num_reviews',
        'branch_dev_enforce_admins',
    ],
];

$pipeline_test_names =
    [
        'has_release' => 'Released',
        'release_after_tools' => 'Released after tools',
        'main_is_release' => 'Main = release',
        'has_json_schema' => 'JSON Schema',
        'has_dsl2_modules_dir' => 'DSL2',
    ] + $base_test_names;
$pipeline_test_descriptions =
    [
        'has_release' => 'Has at least one release',
        'release_after_tools' => 'Last release is after latest tools release (so up to date with template)',
        'main_is_release' => 'Main branch is same commit as the last release',
        'has_json_schema' => 'Has a nextflow_schema.json file (in last release, dev if no release)',
        'has_dsl2_modules_dir' =>
            'Has a modules directory, suggesting that it\'s a DSL2 pipeline (in last release, dev if no release)',
    ] + $base_test_descriptions;
$pipeline_test_descriptions['repo_url'] = 'URL should be set to https://pipelines.tol.sanger.ac.uk/[PIPELINE-NAME]';
$pipeline_test_urls =
    [
        'has_release' => 'https://github.com/sanger-tol/{repo}/releases',
        'release_after_tools' => 'https://github.com/sanger-tol/{repo}/releases/{latest-tag}',
        'main_is_release' => 'https://github.com/sanger-tol/{repo}/compare/{latest-tag}...main',
        'has_json_schema' => 'https://github.com/sanger-tol/{repo}',
        'has_dsl2_modules_dir' => 'https://github.com/sanger-tol/{repo}',
    ] + $base_test_urls;
$pipeline_merge_table_col_headings = $base_merge_table_col_headings;

$core_repo_test_names = $base_test_names;
$core_repo_test_descriptions = $base_test_descriptions;
$core_repo_test_urls = $base_test_urls;
$core_repo_merge_table_col_headings = $base_merge_table_col_headings;
$core_repo_ignore_tests = [
    'branch_dev_exists',
    'branch_template_exists',
    'branch_main_strict_updates',
    'branch_main_required_ci',
    'branch_main_stale_reviews',
    'branch_main_code_owner_reviews',
    'branch_main_required_num_reviews',
    'branch_main_enforce_admins',
    'branch_dev_strict_updates',
    'branch_dev_required_ci',
    'branch_dev_stale_reviews',
    'branch_dev_code_owner_reviews',
    'branch_dev_required_num_reviews',
    'branch_dev_enforce_admins',
    // 'branch_template_restrict_push',
];
foreach ($core_repo_ignore_tests as $key) {
    unset($core_repo_test_names[$key]);
    unset($core_repo_test_descriptions[$key]);
    unset($core_repo_test_urls[$key]);
}

// Get any missing data and run tests / fix problems
$updated_teams = [];
foreach ($core_repos as $idx => $core_repo) {
    $core_repo->test_names = $core_repo_test_names;
    $core_repo->test_descriptions = $core_repo_test_descriptions;
    $core_repo->test_urls = $core_repo_test_urls;
    $core_repo->get_data();
    if ($core_repo->gh_repo->archived) {
        unset($core_repos[$idx]);
        continue;
    }
    $core_repo->run_tests();
    $core_repo->fix_tests();
}
// Tools: Get release info
//$core_repos['tools']->get_release_data();
//$tools_last_release = $core_repos['tools']->gh_release->published_at;
// nf-core tools is not in our sanger-tol organization
$nf_core_tool_repo = new CoreRepoHealth('tools');
$nf_core_tool_repo->get_release_data();
$tools_last_release = $nf_core_tool_repo->gh_release->published_at;

foreach ($pipelines as $idx => $pipeline) {
    $pipeline->test_names = $pipeline_test_names;
    $pipeline->test_descriptions = $pipeline_test_descriptions;
    $pipeline->test_urls = $pipeline_test_urls;
    $pipeline->get_data();
    if ($pipeline->gh_repo->archived) {
        unset($pipelines[$idx]);
        continue;
    }
    $pipeline->run_tests();
    $pipeline->fix_tests();
}

foreach ($updated_teams as $team => $updated) {
    if ($updated) {
        $_GET['action'] = 'refresh';
        get_gh_team_repos($team);
        foreach ($pipelines as $pipeline) {
            $pipeline->test_teams();
        }
        foreach ($core_repos as $core_repo) {
            $core_repo->test_teams();
        }
    }
}

ksort($pipelines);
ksort($core_repos);
?>

<div class="container-fluid main-content">
  <h2>Pipelines</h2>
  <div class="table-responsive">
    <table class="table table-hover table-sm small">
      <thead>
        <tr>
          <th class="small fw-normal text-nowrap">Pipeline Name</th>
          <?php
          $m_names_printed = [];
          $colspan = '';
          foreach ($pipeline_test_names as $key => $name) {
              $description = $pipeline_test_descriptions[$key];
              $print = true;
              foreach ($pipeline_merge_table_col_headings as $m_name => $m_keys) {
                  if (in_array($key, $m_keys)) {
                      if (!in_array($m_name, $m_names_printed)) {
                          $colspan = 'colspan="' . count($m_keys) . '"';
                          $description = $m_name;
                          $name = $m_name;
                          $m_names_printed[] = $m_name;
                      } else {
                          $print = false;
                      }
                  }
              }
              if ($print) {
                  echo '<th ' .
                      $colspan .
                      ' class="small fw-normal text-nowrap" title="' .
                      $description .
                      '" data-bs-toggle="tooltip" data-bs-placement="top">' .
                      $name .
                      '</th>';
              }
          }
          ?>
        </tr>
      </thead>
      <tbody>
        <?php foreach ($pipelines as $pipeline) {
            echo '<tr>';
            echo '<td>' . $pipeline->name . '</td>';
            foreach ($pipeline_test_names as $key => $name) {
                $pipeline->print_table_cell($key);
            }
            echo '</tr>';
        } ?>
      </tbody>
    </table>
  </div>

  <h2>Core repos</h2>
  <div class="table-responsive">
    <table class="table table-hover table-sm small">
      <thead>
        <tr>
          <th class="small text-nowrap">Pipeline Name</th>
          <?php foreach ($core_repo_test_names as $key => $name) {
              echo '<th class="small text-nowrap" title="' .
                  $core_repo_test_descriptions[$key] .
                  '" data-bs-toggle="tooltip" data-bs-placement="top">' .
                  $name .
                  '</th>';
          } ?>
        </tr>
      </thead>
      <tbody>
        <?php foreach ($core_repos as $repo) {
            echo '<tr>';
            echo '<td>' . $repo->name . '</td>';
            foreach ($core_repo_test_names as $key => $name) {
                $repo->print_table_cell($key);
            }
            echo '</tr>';
        } ?>
      </tbody>
    </table>
  </div>

  <h2>Actions</h2>
  <form class="row" action="" method="get">
    <div class="col-4 my-1">
      <select class="form-select repos-select" name="repos">
        <optgroup label="All repositories">
          <option value="all" selected>All pipelines</option>
        </optgroup>
        <optgroup label="Pipelines">
          <?php foreach ($pipelines as $repo) {
              echo '<option>' . $repo->name . '</option>';
          } ?>
        </optgroup>
        <optgroup label="Core Repos">
          <?php foreach ($core_repos as $repo) {
              echo '<option>' . $repo->name . '</option>';
          } ?>
        </optgroup>
      </select>
    </div>
    <button type="submit" name="action" value="refresh" class="btn btn-primary col-2 my-1 ms-2 refresh-btn">Refresh data</button>
    <button type="submit" name="action" value="fix" class="btn btn-info col-2 my-1 ms-2 fix-btn">Fix data</button>
  </form>
  <p><em class="small text-muted">Warning: page will take a minute or two to load. Even when refreshing one repo, some tests will be refreshed for all repos.</em></p>

</div>

<script type="text/javascript">
  $(function() {
    // Disable the buttons to prevent button mashing
    $('.refresh-btn').click(function() {
      $(this).addClass('disabled').html('Refreshing &nbsp; <i class="fas fa-spinner fa-pulse"></i>');
    });
    $('.fix-btn').click(function(e) {
      if (!confirm('This will attempt to change repository settings! Are you sure?')) {
        e.preventDefault();
      } else {
        if ($('.repos-select').val() == 'all') {
          if (!confirm('Seriously - ALL sanger-tol repos. Are you super sure?')) {
            e.preventDefault();
          } else {
            $(this).addClass('disabled').html('Fixing &nbsp; <i class="fas fa-spinner fa-pulse"></i>');
          }
        }
      }
    });

    // Remove all get data from the URL
    if (window.location.href.includes('?')) {
      var url_parts = window.location.href.split('?');
      window.history.replaceState({}, "sanger-tol", url_parts[0]);
    }
  });
</script>



<?php include '../includes/footer.php';
