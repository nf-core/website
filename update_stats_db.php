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

// Get snapshot of key metrics for all repos


$sql = "CREATE TABLE IF NOT EXISTS github_org_users (
                github_org_id     INT AUTO_INCREMENT PRIMARY KEY,
                github_id    INT NOT NULL,
                avatar_url	VARCHAR (400)        DEFAULT NULL,
                github_url	VARCHAR (400)        DEFAULT NULL,
                html_url	VARCHAR (400)        DEFAULT NULL,
                date_added  datetime             DEFAULT current_timestamp
                )";

if (mysqli_query($conn, $sql)) {
    echo "Table created successfully.";
} else {
    echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
}
// Get the current organisation members
// Returns 30 results per page!

$gh_members_url = 'https://api.github.com/orgs/nf-core/members';
$results['gh_org_members'][$updated] = 0;
$first_page = true;
$next_page = false;
while ($first_page || $next_page) {
    // reset loop vars
    $first_page = false;
    // Get GitHub API results
    if ($next_page) {
        $gh_members_url = $next_page;
    }
    $gh_members = json_decode(file_get_contents($gh_members_url, false, $gh_api_opts));
    if (strpos($http_response_header[0], "HTTP/1.1 200") === false) {
        var_dump($http_response_header);
        echo ("\nCould not fetch nf-core members! $gh_members_url");
        continue;
    }


    // Prepare an insert statement
    $sql = "INSERT INTO github_org_users (github_id,avatar_url,github_url,html_url) VALUES (?, ?, ?, ?)";

    if ($stmt = mysqli_prepare($conn, $sql)) {
        // Bind variables to the prepared statement as parameters
        mysqli_stmt_bind_param($stmt, "isss",  $github_id, $avatar_url, $github_url, $html_url,);

        foreach ($gh_members as $idx => $gh_member) {
            $gh_member = get_object_vars($gh_member);
            //check if user already exists
            $check = "SELECT * FROM github_org_users WHERE github_id = '$gh_member[id]'";
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


    $results['gh_org_members'][$updated] += count($gh_members);
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





