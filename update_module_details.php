<?php

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
mysqli_close($conn);
// // Load details of the modules
// $modules_json = json_decode(file_get_contents('public_html/modules.json'));
// $modules = $modules_json->remote_workflows;
// $contribs_try_again = [];

// // Build array of repos to query
// $modules_json_names = [];
// foreach ($modules as $wf) {
//     $modules_json_names[] = $wf->name;
//     if (!isset($results['modules'][$wf->name])) {
//         $results['modules'][$wf->name] = array();
//     }
//     $results['modules'][$wf->name]['num_releases'] = count($wf->releases);
// }
// $ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];
// foreach ($ignored_repos as $name) {
//     if (!isset($results['core_repos'][$name])) {
//         $results['core_repos'][$name] = array();
//     }
// }

// $sql = "CREATE TABLE IF NOT EXISTS github_modules (
//                 id     INT AUTO_INCREMENT PRIMARY KEY,
//                 github_id    INT NOT NULL,
//                 avatar_url	VARCHAR (400)        DEFAULT NULL,
//                 github_url	VARCHAR (400)        DEFAULT NULL,
//                 html_url	VARCHAR (400)        DEFAULT NULL,
//                 date_added  datetime             DEFAULT current_timestamp
//                 )";

// if (mysqli_query($conn, $sql)) {
//     echo "`github_modules` table created successfully.\n";
// } else {
//     echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
// }
// // Get the current organisation members
// // Returns 30 results per page!


// $gh_members = github_query('https://api.github.com/orgs/nf-core/members');
// // Prepare an insert statement
// $sql = "INSERT INTO github_org_members (github_id,avatar_url,github_url,html_url) VALUES (?, ?, ?, ?)";

// if ($stmt = mysqli_prepare($conn, $sql)) {
//     // Bind variables to the prepared statement as parameters
//     mysqli_stmt_bind_param($stmt, "isss",  $github_id, $avatar_url, $github_url, $html_url);

//     foreach ($gh_members as $idx => $gh_member) {
//         //check if user already exists
//         $check = "SELECT * FROM github_org_members WHERE github_id = '$gh_member[id]'";
//         $res = mysqli_query($conn, $check);
//         if ($res->num_rows) {
//             continue;
//         }else{
//             $github_id = $gh_member['id'];
//             $avatar_url = $gh_member['avatar_url'];
//             $github_url = $gh_member['url'];
//             $html_url = $gh_member['html_url'];

//             mysqli_stmt_execute($stmt);
//         }
//     }
// } else {
//     echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
// }

// // github repos
// $sql = "CREATE TABLE IF NOT EXISTS github_repos(
//                 id                  INT AUTO_INCREMENT PRIMARY KEY,
//                 github_id           INT NOT NULL,
//                 name	            VARCHAR (400)    DEFAULT NULL,
//                 default_branch  	VARCHAR (400)    DEFAULT NULL,
//                 size                INT              DEFAULT NULL,
//                 stargazers_count    INT              DEFAULT NULL,
//                 forks_count         INT              DEFAULT NULL,
//                 open_issues_count   INT              DEFAULT NULL,
//                 date_created        datetime         DEFAULT NULL,
//                 date_pushed         datetime         DEFAULT NULL
//                 )";

// if (mysqli_query($conn, $sql)) {
//     echo "'github_repos' table created successfully.\n";
// } else {
//     echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
// }
// $gh_repos = github_query('https://api.github.com/orgs/nf-core/repos');

// $sql = "INSERT INTO github_repos (github_id, name, default_branch, size, stargazers_count, forks_count, open_issues_count, date_created, date_pushed) VALUES (?, ?, ?, ?,?, ?, ?, ?,?)";

// if ($stmt = mysqli_prepare($conn, $sql)) {
//     // Bind variables to the prepared statement as parameters
//     mysqli_stmt_bind_param($stmt, "issiiiiss",  $github_id , $name, $default_branch, $size, $stargazers_count, $forks_count, $open_issues_count, $date_created, $date_pushed);
//     $ignored_repos = parse_ini_file("ignored_repos.ini")['repos'];
    
//     foreach ($gh_repos as $idx => $gh_repo) {
//         if(in_array($gh_repo['name'],$ignored_repos)){
//             continue;
//         }
//         //check if user already exists
//         $check = "SELECT * FROM github_repos WHERE github_id = '$gh_repo[id]'";
//         $res = mysqli_query($conn, $check);
//         if ($res->num_rows) {
//             continue;
//         } else {
//             $github_id  = $gh_repo['id'];
//             $name = $gh_repo['name'];
//             $default_branch = $gh_repo['default_branch'];
//             $size = $gh_repo['size'];
//             $stargazers_count = $gh_repo['stargazers_count'];
//             $forks_count = $gh_repo['forks_count'];
//             $open_issues_count = $gh_repo['open_issues_count'];
//             $date_created = date_format(new DateTime($gh_repo['created_at']), 'Y-m-d H:i:s');
//             $date_pushed = date_format(new DateTime($gh_repo['pushed_at']), 'Y-m-d H:i:s');

//             $inserted = mysqli_stmt_execute($stmt);
//             if(!$inserted){
//                 echo "ERROR: " . mysqli_stmt_error($stmt);
//             }
//         }
//     }
// } else {
//     echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
// }
// $check = "SELECT id, name FROM github_repos";

// if ($result = mysqli_query($conn, $check)) {
//     if (mysqli_num_rows($result) > 0) {
//         $repos = mysqli_fetch_all($result);
        
//     } else {
//         echo "ERROR: No entires in 'github_repos'";
//     }
// } else {
//     echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
// }
// // github commits
// $sql = "CREATE TABLE IF NOT EXISTS github_commits(
//                 id                  INT AUTO_INCREMENT PRIMARY KEY,
//                 sha                 VARCHAR (400)    NOT NULL,
//                 author_id           INT,
//                 repo_id             INT,
//                 message             VARCHAR (800)    DEFAULT NULL,
//                 date                datetime         DEFAULT NULL
//                 )";

// if (mysqli_query($conn, $sql)) {
//     echo "'github_commits' table created successfully.\n";
// } else {
//     echo "ERROR: Could not able to execute $sql. " . mysqli_error($conn);
// }

// foreach ($repos as $key => $repo) {
//     $gh_commits = github_query('https://api.github.com/repos/nf-core/'.$repo[1].'/commits');
//     print_r($repo[1] . "\n");
//     $sql = "INSERT INTO github_commits (sha, author_id, repo_id, message, date) VALUES (?, ?, ?, ?, ?)";

//     if ($stmt = mysqli_prepare($conn, $sql)) {
//         // Bind variables to the prepared statement as parameters
//         mysqli_stmt_bind_param($stmt, "siiss",$sha, $author_id, $repo_id, $message, $date);
        
//         foreach ($gh_commits as $idx => $gh_commit) {
//             //check if commit already exists
//             $check = "SELECT * FROM github_commits WHERE sha = '$gh_commit[sha]'";
//             $res = mysqli_query($conn, $check);
//             if ($res->num_rows) {
//                 continue;
//             } else {
//                 $message_trunc = $gh_commit['commit']['message'];
//                 if (strlen($message_trunc) > 800){
//                     $message_trunc = substr($message_trunc, 0, 796) . "...";
//                 } 
//                 if (!$gh_commit['author']){
//                     #try and guess username
//                     $author_guess = github_query("https://api.github.com/search/users?q=". urlencode($gh_commit['commit']['author']['name']));
//                     print_r($author_guess);
//                     if($author_guess['total_count']>0){
//                         $author_id = $author_guess['items'][0]['id'];
//                     }  
//                 } else {
//                     $author_id = $gh_commit['author']['id'];
//                 }

//                 $sha = $gh_commit['sha'];
                
//                 $repo_id = $repo[0];
//                 $message = $message_trunc;
//                 $date= date_format(new DateTime($gh_commit['commit']['author']['date']), 'Y-m-d H:i:s');

//                 $inserted = mysqli_stmt_execute($stmt);
//                 if (!$inserted) {
//                     echo "ERROR: " . mysqli_stmt_error($stmt);
//                 }
//             }
//         }
//     } else {
//         echo "ERROR: Could not prepare query: $sql. " . mysqli_error($conn);
//     }
// }
