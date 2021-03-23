<?php

// Allow PHP fopen to work with remote links
ini_set("allow_url_fopen", 1);

$updated = time();

// Get auth secrets
$config = parse_ini_file("config.ini");
$gh_auth = base64_encode($config['github_username'] . ':' . $config['github_access_token']);

$config = parse_ini_file("config.ini");
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

        array_push($res, ...$tmp_results);

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

// Load details of the pipelines
$pipelines_json = json_decode(file_get_contents('public_html/pipelines.json'));
$pipelines = $pipelines_json->remote_workflows;
$contribs_try_again = [];

// Build array of repos to query
$pipelines_json_names = [];
foreach ($pipelines as $wf) {
    $pipelines_json_names[] = $wf->name;
    if (!isset($results['pipelines'][$wf->name])) {
        $results['pipelines'][$wf->name] = array();
    }
    $results['pipelines'][$wf->name]['num_releases'] = count($wf->releases);
}
$ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];
foreach ($ignored_repos as $name) {
    if (!isset($results['core_repos'][$name])) {
        $results['core_repos'][$name] = array();
    }
}

$sql = "CREATE TABLE IF NOT EXISTS github_org_members (
                id     INT AUTO_INCREMENT PRIMARY KEY,
                github_id    INT NOT NULL,
                avatar_url	VARCHAR (400)        DEFAULT NULL,
                github_url	VARCHAR (400)        DEFAULT NULL,
                html_url	VARCHAR (400)        DEFAULT NULL,
                date_added  datetime             DEFAULT current_timestamp
                )";

if (mysqli_query($conn, $sql)) {
    echo "`github_org_members` table created successfully.\n";
} else {
    echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
}
// Get the current organisation members
// Returns 30 results per page!


$gh_members = github_query('https://api.github.com/orgs/nf-core/members');
// Prepare an insert statement
$sql = "INSERT INTO github_org_members (github_id,avatar_url,github_url,html_url) VALUES (?, ?, ?, ?)";

if ($stmt = mysqli_prepare($conn, $sql)) {
    // Bind variables to the prepared statement as parameters
    mysqli_stmt_bind_param($stmt, "isss",  $github_id, $avatar_url, $github_url, $html_url,);

    foreach ($gh_members as $idx => $gh_member) {
        //check if user already exists
        $check = "SELECT * FROM github_org_members WHERE github_id = '$gh_member[id]'";
        $res = mysqli_query($conn, $check);
        if ($res->num_rows) {
            continue;
        }else{
            $github_id = $gh_member['id'];
            $avatar_url = $gh_member['avatar_url'];
            $github_url = $gh_member['url'];
            $html_url = $gh_member['html_url'];

            mysqli_stmt_execute($stmt);
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}

//
//
// SLACK API calls
//
//
function slack_query($url){
    global $config;

    $slack_api_opts = stream_context_create([
        'http' => [
            'method' => 'GET',
            'header' => [
                'Authorization: Bearer ' . $config['slack_access_token'],

            ]
        ]
    ]);
    $res = json_decode(file_get_contents($url, false, $slack_api_opts), true);
    $next_page = $res['response_metadata']['next_cursor'];

    while ($next_page != "") {
        $url .= '&cursor=' . urlencode($next_page);
        $tmp_results = json_decode(file_get_contents($url, false, $slack_api_opts), true);
        array_push($res['members'], ...$tmp_results['members']);
        $next_page = $tmp_results['response_metadata']['next_cursor'];
    }
    return $res;
}

// Slack users
$sql = "CREATE TABLE IF NOT EXISTS slack_users (
                id          INT AUTO_INCREMENT PRIMARY KEY,
                slack_id    VARCHAR(20)          NOT NULL,
                user_name	VARCHAR (400)        DEFAULT NULL,
                first_name	VARCHAR (400)        DEFAULT NULL,
                last_name	VARCHAR (400)        DEFAULT NULL,
                avatar_url	VARCHAR (400)        DEFAULT NULL,
                time_zone	VARCHAR (400)        DEFAULT NULL,
                date_added  datetime             DEFAULT current_timestamp
                )";

if (mysqli_query($conn, $sql)) {
    echo "`slack_users` Table created successfully. \n";
} else {
    echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
}
$slack_users = slack_query('https://slack.com/api/users.list?pretty=1');
$sql = "INSERT INTO slack_users(slack_id,user_name,first_name,last_name,avatar_url,time_zone) VALUES (?, ?, ?, ?, ?, ?)";

if ($stmt = mysqli_prepare($conn, $sql)) {
// Bind variables to the prepared statement as parameters
mysqli_stmt_bind_param($stmt, "ssssss", $slack_id, $user_name, $first_name, $last_name, $avatar_url, $time_zone);

    foreach ($slack_users['members'] as $idx => $slack_user) {
        //check if user already exists
        $check = "SELECT * FROM slack_users WHERE slack_id = '$slack_user[id]'";
        $res = mysqli_query($conn, $check);
        if ($res->num_rows) {
            continue; //user exists already
        } else {
            $slack_id = $slack_user['id'];
            $user_name = $slack_user['name'];
            $first_name = isset($slack_user['profile']['first_name']) ? $slack_user['profile']['first_name'] : NULL;
            $last_name = isset($slack_user['profile']['last_name']) ? $slack_user['profile']['last_name'] : NULL;
            $avatar_url = isset($slack_user['profile']['image_512']) ? $slack_user['profile']['image_512'] : NULL;
            $time_zone = isset($slack_user['tz']) ? $slack_user['tz'] : NULL;
            mysqli_stmt_execute($stmt);
        }
    }
} else {
    echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
}
