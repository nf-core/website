<?php
echo "\n\nUpdating pipeline details - " . date("Y-m-d h:i:s") . "\n";

// Allow PHP fopen to work with remote links
ini_set("allow_url_fopen", 1);
// Load yaml parser
require "vendor/autoload.php";
use Symfony\Component\Yaml\Yaml;

// Get auth secrets
$config = parse_ini_file("config.ini");
$gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

if ($conn === false) {
    die("ERROR: Could not connect. " . mysqli_connect_error());
}

//
//
// GitHub API calls
//
//

function github_query($gh_query_url)
{
    global $config;
    $gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

    // HTTP header to use on GitHub API GET requests
    $gh_api_opts = stream_context_create([
        'http' => [
            'method' => 'GET',
            'header' => [
                'User-Agent: PHP',
                "Authorization: Basic $gh_auth"
            ]
        ]
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
        if (strpos($http_response_header[0], "HTTP/1.1 200") === false) {
            var_dump($http_response_header);
            echo ("\nCould not fetch $gh_query_url");
            continue;
        }
        if(substr($gh_query_url,0,29)== "https://api.github.com/repos/" ){
            $res = $tmp_results;
        } else{
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

$gh_modules = github_query('https://api.github.com/repos/nf-core/modules/git/trees/master?recursive=1');
$modules = [];
foreach($gh_modules['tree'] as $f){
    if(substr($f['path'],-8)=='meta.yml' && substr($f['path'],0,8)=='modules/'){
        $meta = github_query($f['url']);
        $meta_content = base64_decode($meta['content']);
        $meta_content = Yaml::parse($meta_content);
        $meta_content['keywords'] = is_array($meta_content['keywords']) ? implode(';', $meta_content['keywords']) : $meta_content['keywords'];
        $meta_content['authors'] = is_array($meta_content['authors']) ? implode(';', $meta_content['authors']) : $meta_content['authors'];
        $meta_content['sha'] = $meta['sha'];
        $meta_content['api_url'] = $meta['url'];
        $meta_content['github_path'] = $f['path'];
        $modules[] = $meta_content;
    }
    
}
// Drop existing table if query was successful
if (count($modules)>1){
    $sql = "DROP TABLE IF EXISTS nfcore_modules";
    if (!mysqli_query($conn, $sql)) {
        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
    }
    $sql = "CREATE TABLE nfcore_modules (
                id                          INT AUTO_INCREMENT PRIMARY KEY,
                github_sha  VARCHAR (400)   NOT NULL,
                github_path VARCHAR (400)   NOT NULL,
                api_url     VARCHAR (400)   NOT NULL,
                name	    VARCHAR (400)   NOT NULL,
                description	VARCHAR (4000)  DEFAULT NULL,
                keywords	VARCHAR (2000)   DEFAULT NULL,
                tools	    JSON            NOT NULL,
                input	    JSON            NOT NULL,
                output	    JSON            NOT NULL,
                authors	    VARCHAR (2000)   NOT NULL,
                date_added  datetime        DEFAULT current_timestamp
                )";
    if (mysqli_query($conn, $sql)) {
        echo "`nfcore_modules` table created successfully.\n";
    } else {
        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
    }
}
// Prepare an insert statement
$sql = "INSERT INTO nfcore_modules (github_sha,github_path,api_url,name,description,keywords,tools,input,output,authors) VALUES (?,?,?,?,?,?,?,?,?,?)";

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param($stmt, "ssssssssss", $github_sha, $github_path, $api_url, $name, $description, $keywords, $tools, $input, $output, $authors);

    foreach ($modules as $idx => $module) {
        // check if user already exists
        $check = "SELECT * FROM nfcore_modules WHERE name = '$module[name]'";
        $res = mysqli_query($conn, $check);
        if ($res->num_rows) {
            continue;
        }else{
            $github_sha = $module['sha'];
            $github_path = $module['github_path'];
            $api_url = $module['api_url'];
            $name = $module['name'];
            $description = $module['description'];
            $keywords = $module['keywords'];
            $tools = json_encode($module['tools']);
            $input = json_encode($module['input']);
            $output = json_encode($module['output']);
            $authors =  $module['authors'];

            mysqli_stmt_execute($stmt);
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}

//
//  nf-core pipelines table
//

$gh_pipelines = github_query('https://api.github.com/orgs/nf-core/repos?per_page=100');
$ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];

// Drop existing table if query was successful
if (count($gh_pipelines) > 1) {
    $sql = "DROP TABLE IF EXISTS nfcore_pipelines";
    if (!mysqli_query($conn, $sql)) {
        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
    }
    $sql = "CREATE TABLE nfcore_pipelines (
                id          INT             AUTO_INCREMENT PRIMARY KEY,
                github_id  VARCHAR (400)   NOT NULL,
                html_url     VARCHAR (400)   NOT NULL,
                name	    VARCHAR (400)   NOT NULL,
                description	VARCHAR (4000)  DEFAULT NULL,
                gh_created_at datetime       NOT NULL,
                gh_updated_at datetime       NOT NULL,
                gh_pushed_at datetime       NOT NULL,
                stargazers_count INT         NOT NULL,
                watchers_count INT           NOT NULL,
                forks_count INT             NOT NULL,
                open_issues_count INT       NOT NULL,
                topics       VARCHAR (4000)  DEFAULT NULL,
                watchers     INT             NOT NULL,
                default_branch VARCHAR (400) NOT NULL,
                pipeline_type VARCHAR (400)  NOT NULL,
                date_added  datetime        DEFAULT current_timestamp
                )";
    if (mysqli_query($conn, $sql)) {
        echo "`nfcore_pipelines` table created successfully.\n";
    } else {
        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
    }
}
// Prepare an insert statement
$sql = "INSERT INTO nfcore_pipelines (github_id,html_url,name,description,gh_created_at,gh_updated_at,gh_pushed_at,stargazers_count,watchers_count,forks_count,open_issues_count,topics,watchers,default_branch,pipeline_type) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param($stmt, "sssssssiiiisiss",$github_id,$html_url,$name,$description,$gh_created_at,$gh_updated_at,$gh_pushed_at,$stargazers_count,$watchers_count,$forks_count,$open_issues_count,$topics,$watchers,$default_branch,$pipeline_type);
    foreach ($gh_pipelines as $pipeline) {
        // check if user already exists
        
        
            $github_id = $pipeline['id'];
            $html_url = $pipeline['html_url'];
            $name = $pipeline['name'];
            $description = $pipeline['description'];
            $gh_created_at = date('Y-m-d H:i:s',$pipeline['created_at']);
            $gh_updated_at = date('Y-m-d H:i:s',$pipeline['updated_at']);
            $gh_pushed_at = date('Y-m-d H:i:s',$pipeline['pushed_at']);
            $stargazers_count = $pipeline['stargazers_count'];
            $watchers_count = $pipeline['watchers_count'];
            $forks_count = $pipeline['forks_count'];
            $open_issues_count = $pipeline['open_issues_count'];
            $topics = is_array($pipeline['topics']) ? implode(';', $pipeline['topics']) : $pipeline['topics'];
            $watchers = $pipeline['watchers'];
            $default_branch = $pipeline['default_branch'];

            if (in_array($pipeline['name'], $ignored_repos)) {
                $pipeline_type = 'core_repos';
            } else {
                $pipeline_type = 'pipelines';
            }
            mysqli_stmt_execute($stmt);
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
        echo "Oops! Something went wrong. Please try again later.";
    }
}
// Drop existing table if query was successful
if (count($pipelines) > 1) {
    $sql = "DROP TABLE IF EXISTS pipelines_modules";
    if (!mysqli_query($conn, $sql)) {
        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
    }
    $sql = "CREATE TABLE nfcore_pipelines (
                id          INT             AUTO_INCREMENT PRIMARY KEY,
                github_id  VARCHAR (400)   NOT NULL,
                html_url     VARCHAR (400)   NOT NULL,
                name	    VARCHAR (400)   NOT NULL,
                description	VARCHAR (4000)  DEFAULT NULL,
                gh_created_at datetime       NOT NULL,
                gh_updated_at datetime       NOT NULL,
                gh_pushed_at datetime       NOT NULL,
                stargazers_count INT         NOT NULL,
                watchers_count INT           NOT NULL,
                forks_count INT             NOT NULL,
                open_issues_count INT       NOT NULL,
                topics       VARCHAR (4000)  DEFAULT NULL,
                watchers     INT             NOT NULL,
                default_branch VARCHAR (400) NOT NULL,
                pipeline_type VARCHAR (400)  NOT NULL,
                date_added  datetime        DEFAULT current_timestamp
                )";
    if (mysqli_query($conn, $sql)) {
        echo "`nfcore_pipelines` table created successfully.\n";
    } else {
        echo "ERROR: Could not execute $sql. " . mysqli_error($conn);
    }
}
// Prepare an insert statement
$sql = "INSERT INTO nfcore_pipelines (github_id,html_url,name,description,gh_created_at,gh_updated_at,gh_pushed_at,stargazers_count,watchers_count,forks_count,open_issues_count,topics,watchers,default_branch,pipeline_type) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param($stmt, "sssssssiiiisiss", $github_id, $html_url, $name, $description, $gh_created_at, $gh_updated_at, $gh_pushed_at, $stargazers_count, $watchers_count, $forks_count, $open_issues_count, $topics, $watchers, $default_branch, $pipeline_type);
    foreach ($gh_pipelines as $pipeline) {
        // check if user already exists


        $github_id = $pipeline['id'];
        $html_url = $pipeline['html_url'];
        $name = $pipeline['name'];
        $description = $pipeline['description'];
        $gh_created_at = date('Y-m-d H:i:s', $pipeline['created_at']);
        $gh_updated_at = date('Y-m-d H:i:s', $pipeline['updated_at']);
        $gh_pushed_at = date('Y-m-d H:i:s', $pipeline['pushed_at']);
        $stargazers_count = $pipeline['stargazers_count'];
        $watchers_count = $pipeline['watchers_count'];
        $forks_count = $pipeline['forks_count'];
        $open_issues_count = $pipeline['open_issues_count'];
        $topics = is_array($pipeline['topics']) ? implode(';', $pipeline['topics']) : $pipeline['topics'];
        $watchers = $pipeline['watchers'];
        $default_branch = $pipeline['default_branch'];

        if (in_array($pipeline['name'], $ignored_repos)) {
            $pipeline_type = 'core_repos';
        } else {
            $pipeline_type = 'pipelines';
        }
        mysqli_stmt_execute($stmt);
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}

mysqli_close($conn);
echo ("\nupdate_pipeline_details done " . date("Y-m-d h:i:s") . "\n\n");
