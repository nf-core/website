<?php
echo "\n\nUpdating module details - " . date('Y-m-d h:i:s') . "\n";

// Allow PHP fopen to work with remote links
ini_set('allow_url_fopen', 1);
// Load yaml parser
require 'vendor/autoload.php';
use Symfony\Component\Yaml\Yaml;

// Get auth secrets
$config = parse_ini_file('config.ini');
$gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

if ($conn === false) {
    die('ERROR: Could not connect. ' . mysqli_connect_error());
}

//
//
// GitHub API calls
//
//

function github_query($gh_query_url) {
    global $config;
    $gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

    // HTTP header to use on GitHub API GET requests
    $gh_api_opts = stream_context_create([
        'http' => [
            'method' => 'GET',
            'header' => ['User-Agent: PHP', "Authorization: Basic $gh_auth"],
        ],
    ]);

    $first_page = true;
    $next_page = false;
    $res = [];
    while ($first_page || $next_page) {
        // reset loop vars
        $first_page = false;
        // Get GitHub API results
        if ($next_page) {
            $gh_query_url = $next_page;
        }
        $tmp_results = json_decode(file_get_contents($gh_query_url, false, $gh_api_opts), true);
        if (strpos($http_response_header[0], 'HTTP/1.1 200') === false) {
            var_dump($http_response_header);
            echo "\nCould not fetch $gh_query_url";
            continue;
        }
        if (substr($gh_query_url, 0, 29) == 'https://api.github.com/repos/') {
            $res = $tmp_results;
        } else {
            array_push($res, ...$tmp_results);
        }
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
    return $res;
}

//
//  nf-core modules table
//

$gh_modules = github_query('https://api.github.com/repos/sanger-tol/nf-core-modules/git/trees/main?recursive=1');
$modules = [];
foreach ($gh_modules['tree'] as $f) {
    if (substr($f['path'], -8) == 'meta.yml' && substr($f['path'], 0, 8) == 'modules/') {
        $meta = github_query($f['url']);
        $meta_content = base64_decode($meta['content']);
        $meta_content = Yaml::parse($meta_content);
        $meta_content['keywords'] = is_array($meta_content['keywords'])
            ? implode(';', $meta_content['keywords'])
            : $meta_content['keywords'];
        $meta_content['authors'] = is_array($meta_content['authors'])
            ? implode(';', $meta_content['authors'])
            : $meta_content['authors'];
        $meta_content['sha'] = $meta['sha'];
        $meta_content['api_url'] = $meta['url'];
        $meta_content['github_path'] = $f['path'];
        $modules[] = $meta_content;
    }
}

$sql = "CREATE TABLE IF NOT EXISTS  nfcore_modules (
            id           INT             AUTO_INCREMENT PRIMARY KEY,
            github_sha   VARCHAR (400)   NOT NULL,
            github_path  VARCHAR (400)   NOT NULL,
            api_url      VARCHAR (400)   NOT NULL,
            name	     VARCHAR (400)   NOT NULL,
            description	 VARCHAR (4000)  DEFAULT NULL,
            keywords	 VARCHAR (2000)  DEFAULT NULL,
            tools	     JSON            NOT NULL,
            input	     JSON            NOT NULL,
            output	     JSON            NOT NULL,
            authors	     VARCHAR (2000)  NOT NULL,
            date_added   TIMESTAMP       DEFAULT CURRENT_TIMESTAMP,
            date_updated TIMESTAMP       DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
            )";
if (mysqli_query($conn, $sql)) {
    echo "`nfcore_modules` table created successfully.\n";
} else {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}
// Prepare an insert statement
$sql =
    'INSERT INTO nfcore_modules (github_sha,github_path,api_url,name,description,keywords,tools,input,output,authors) VALUES (?,?,?,?,?,?,?,?,?,?)';

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param(
        $stmt,
        'ssssssssss',
        $github_sha,
        $github_path,
        $api_url,
        $name,
        $description,
        $keywords,
        $tools,
        $input,
        $output,
        $authors,
    );

    foreach ($modules as $idx => $module) {
        // check if module already exists
        $check = "SELECT * FROM nfcore_modules WHERE name = '$module[name]'";
        $res = mysqli_query($conn, $check);
        $github_sha = $module['sha'];
        $github_path = $module['github_path'];
        $api_url = $module['api_url'];
        $name = $module['name'];
        $description = $module['description'];
        $keywords = $module['keywords'];
        $tools = json_encode($module['tools']);
        $input = json_encode($module['input']);
        $output = json_encode($module['output']);
        $authors = $module['authors'];
        if ($res->num_rows) {
            // check if entry is different
            $row = $res->fetch_assoc();

            // update existing module
            $update = "UPDATE nfcore_modules SET github_sha = ?, github_path = ?, api_url = ?, description = ?, keywords = ?, tools = ?, input = ?, output = ?, authors = ? WHERE name = '$name'";
            if ($update_stmt = mysqli_prepare($conn, $update)) {
                // Bind variables to the prepared statement as parameters
                mysqli_stmt_bind_param(
                    $update_stmt,
                    'sssssssss',
                    $github_sha,
                    $github_path,
                    $api_url,
                    $description,
                    $keywords,
                    $tools,
                    $input,
                    $output,
                    $authors,
                );
            } else {
                echo "ERROR: Could not prepare query: $update. " . mysqli_error($conn);
            }
            if (!mysqli_stmt_execute($update_stmt)) {
                echo "ERROR: Could not execute $update. " . mysqli_error($conn);
            }
        } else {
            if (!mysqli_stmt_execute($stmt)) {
                echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
            }
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}

//
//  nf-core pipelines table
//

$gh_pipelines = github_query('https://api.github.com/orgs/sanger-tol/repos?per_page=100');
$ignored_repos = parse_ini_file('ignored_repos.ini')['repos'];

// Drop existing table if query was successful
$sql = "CREATE TABLE IF NOT EXISTS nfcore_pipelines (
            id                INT             AUTO_INCREMENT PRIMARY KEY,
            github_id         VARCHAR (400)   NOT NULL,
            html_url          VARCHAR (400)   NOT NULL,
            name	          VARCHAR (400)   NOT NULL,
            description	      VARCHAR (4000)  DEFAULT NULL,
            gh_created_at     datetime        NOT NULL,
            gh_updated_at     datetime        NOT NULL,
            gh_pushed_at      datetime        NOT NULL,
            stargazers_count  INT             NOT NULL,
            watchers_count    INT             NOT NULL,
            forks_count       INT             NOT NULL,
            open_issues_count INT             NOT NULL,
            open_pr_count     INT             NOT NULL,
            topics            VARCHAR (4000)  DEFAULT NULL,
            default_branch    VARCHAR (400)   NOT NULL,
            pipeline_type     VARCHAR (400)   NOT NULL,
            archived          BOOLEAN         NOT NULL,
            last_release_date datetime        DEFAULT NULL,
            date_added        datetime        DEFAULT current_timestamp
            )";
if (mysqli_query($conn, $sql)) {
    echo "`nfcore_pipelines` table created successfully.\n";
} else {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}
// Prepare an insert statement
$sql =
    'INSERT INTO nfcore_pipelines (github_id,html_url,name,description,gh_created_at,gh_updated_at,gh_pushed_at,stargazers_count,watchers_count,forks_count,open_issues_count,open_pr_count,topics,default_branch,pipeline_type,archived,last_release_date) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)';

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param(
        $stmt,
        'sssssssiiiiisssis',
        $github_id,
        $html_url,
        $name,
        $description,
        $gh_created_at,
        $gh_updated_at,
        $gh_pushed_at,
        $stargazers_count,
        $watchers_count,
        $forks_count,
        $open_issues_count,
        $open_pr_count,
        $topics,
        $default_branch,
        $pipeline_type,
        $archived,
        $last_release_date,
    );
    foreach ($gh_pipelines as $pipeline) {
        // check where entries need to be updated and update them

        $github_id = $pipeline['id'];
        $html_url = $pipeline['html_url'];
        $name = $pipeline['name'];
        $description = $pipeline['description'];
        $gh_created_at = date('Y-m-d H:i:s', strtotime($pipeline['created_at']));
        $gh_updated_at = date('Y-m-d H:i:s', strtotime($pipeline['updated_at']));
        $gh_pushed_at = date('Y-m-d H:i:s', strtotime($pipeline['pushed_at']));
        $stargazers_count = $pipeline['stargazers_count'];
        $watchers_count = count(
            github_query('https://api.github.com/repos/sanger-tol/' . $pipeline['name'] . '/watchers'),
        );
        $forks_count = $pipeline['forks_count'];
        $open_issues_count = $pipeline['open_issues_count'];
        $open_pr_count = count(github_query(str_replace('{/number}', '', $pipeline['pulls_url'])));
        $topics = is_array($pipeline['topics']) ? implode(';', $pipeline['topics']) : $pipeline['topics'];
        $default_branch = $pipeline['default_branch'];
        $archived = $pipeline['archived'];
        $last_release_date = github_query(str_replace('{/id}', '', $pipeline['releases_url']) . '?per_page=1')[0][
            'published_at'
        ];
        $last_release_date = is_null($last_release_date) ? null : date('Y-m-d H:i:s', strtotime($last_release_date));
        if (in_array($pipeline['name'], $ignored_repos)) {
            $pipeline_type = 'core_repos';
        } else {
            $pipeline_type = 'pipelines';
        }
        $check = "SELECT * FROM nfcore_pipelines WHERE name = '" . $pipeline['name'] . "'";
        $res = mysqli_query($conn, $check);
        if ($res->num_rows > 0) {
            // test if the entry needs to be updated
            $row = $res->fetch_assoc();
            $update = false;
            $update = $update || $row['github_id'] != $github_id;
            $update = $update || $row['html_url'] != $html_url;
            $update = $update || $row['description'] != $description;
            $update = $update || $row['gh_created_at'] != $gh_created_at;
            $update = $update || $row['gh_updated_at'] != $gh_updated_at;
            $update = $update || $row['gh_pushed_at'] != $gh_pushed_at;
            $update = $update || $row['stargazers_count'] != $stargazers_count;
            $update = $update || $row['watchers_count'] != $watchers_count;
            $update = $update || $row['forks_count'] != $forks_count;
            $update = $update || $row['open_issues_count'] != $open_issues_count;
            $update = $update || $row['topics'] != $topics;
            $update = $update || $row['default_branch'] != $default_branch;
            $update = $update || $row['pipeline_type'] != $pipeline_type;
            $update = $update || $row['archived'] != $archived;
            $update = $update || $row['last_release_date'] != $last_release_date;
            if ($update) {
                $update = 'UPDATE nfcore_pipelines SET ';
                $update .= "github_id =  '$github_id',";
                $update .= "html_url =  '$html_url',";
                $update .= "description =  '$description',";
                $update .= "gh_created_at =  '$gh_created_at',";
                $update .= "gh_updated_at =  '$gh_updated_at',";
                $update .= "gh_pushed_at =  '$gh_pushed_at',";
                $update .= "stargazers_count =  '$stargazers_count',";
                $update .= "watchers_count =  '$watchers_count',";
                $update .= "forks_count =  '$forks_count',";
                $update .= "open_issues_count =  '$open_issues_count',";
                $update .= "topics =  '$topics',";
                $update .= "default_branch =  '$default_branch',";
                $update .= "pipeline_type =  '$pipeline_type',";
                $update .= "archived =  '$archived',";
                $update .= "last_release_date =  '$last_release_date'";
                $update .= " WHERE name = '" . $pipeline['name'] . "'";

                if (mysqli_query($conn, $update)) {
                    echo "Updated $pipeline[name]\n";
                } else {
                    mysqli_error($conn);
                }
            }
        } else {
            if (!mysqli_stmt_execute($stmt)) {
                echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
            }
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}

//
//  pipelines modules link table
//

$sql = "SELECT * FROM nfcore_pipelines WHERE pipeline_type = 'pipelines' ORDER BY LOWER(name) ";
$pipelines = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        $pipelines = mysqli_fetch_all($result, MYSQLI_ASSOC);
        // Free result set
        mysqli_free_result($result);
    } else {
        echo 'No pipelines found.';
    }
}

// create pipelines_modules table with foreign keys to pipelines and modules
$sql = "CREATE TABLE IF NOT EXISTS  pipelines_modules (
                id          INT             AUTO_INCREMENT PRIMARY KEY,
                pipeline_id INT             NOT NULL,
                module_id   INT             NOT NULL,
                FOREIGN KEY (pipeline_id)   REFERENCES nfcore_pipelines(id),
                FOREIGN KEY (module_id)     REFERENCES nfcore_modules(id)
                )";
if (mysqli_query($conn, $sql)) {
    echo "`pipelines_modules` table created successfully.\n";
} else {
    echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
}
// Prepare an insert statement
$sql = 'INSERT INTO pipelines_modules (pipeline_id,module_id) VALUES (?,?)';

foreach ($pipelines as $pipeline) {
    $modules_json = github_query(
        'https://api.github.com/repos/sanger-tol/' . $pipeline['name'] . '/contents/modules.json',
    );
    $modules_json = json_decode(base64_decode($modules_json['content']), true);

    /*
    How our modules.json file being configed, use nf-core or sanger-tol?
    Not well formated and has no startard format
    */

    $modules = $modules_json['repos']['nf-core/modules'];
 

    // catch repos with no modules.json
    if ($modules == null) {
        continue;
    }
    foreach ($modules as $name => $content) {
        $stmt = mysqli_prepare($conn, $sql);
        // Bind variables to the prepared statement as parameters
        mysqli_stmt_bind_param($stmt, 'ii', $pipeline_id, $module_id);

        $name = str_replace('/', '_', $name);
        // pepare a select statment for nfcore_modules based on name
        $get_module = "SELECT * FROM nfcore_modules WHERE name = '$name'";

        if ($result = mysqli_query($conn, $get_module)) {
            if (mysqli_num_rows($result) > 0) {
                $repo_module = mysqli_fetch_all($result, MYSQLI_ASSOC);
                $pipeline_id = $pipeline['id'];
                $module_id = $repo_module[0]['id'];

                # need to check before inserting
                $check_sql = "SELECT * FROM pipelines_modules WHERE pipeline_id = '$pipeline_id' AND module_id  = '$module_id'";
                $res = mysqli_query($conn, $check_sql);
                if ($res->num_rows == 0) {                
                    mysqli_stmt_execute($stmt);
                    print_r($stmt->error);
                }
                // Free result set
                mysqli_free_result($result);
            } else {
                echo "No modules with the name $name found for pipeline ". $pipeline['name']."\n";
            }
        }
    }
}

mysqli_close($conn);
echo "\nupdate_module_details done " . date('Y-m-d h:i:s') . "\n\n";
