<?php
// Outgoing HTTP(S) requests through the curl extension, so they keep working on
// hosts where allow_url_fopen is disabled (e.g. ea-php84 on the oldsite server).
//
// fetch_url($url, $context) is a drop-in for file_get_contents($url, false, $context):
// it accepts the same stream context (method, header, content, timeout, ignore_errors)
// and returns the response body, or false on a network error or an HTTP error status.
// fetch_url_headers() returns the response headers of the last request, status line
// first ("HTTP/1.1 200 OK"), like the old $http_response_header.

function fetch_url_headers(?array $set = null): array {
    static $headers = [];
    if ($set !== null) {
        $headers = $set;
    }
    return $headers;
}

function fetch_url(string $url, $context = null) {
    // Without curl, fall back to the stream wrapper (needs allow_url_fopen)
    if (!extension_loaded('curl')) {
        $body = $context ? file_get_contents($url, false, $context) : file_get_contents($url);
        fetch_url_headers(http_get_last_response_headers() ?? []);
        return $body;
    }

    $http = $context ? stream_context_get_options($context)['http'] ?? [] : [];
    $method = strtoupper($http['method'] ?? 'GET');
    $request_headers = $http['header'] ?? [];
    if (is_string($request_headers)) {
        $request_headers = preg_split('/\r\n|\n/', trim($request_headers));
    }

    $response_headers = [];
    $ch = curl_init($url);
    curl_setopt_array($ch, [
        CURLOPT_RETURNTRANSFER => true,
        CURLOPT_FOLLOWLOCATION => true,
        CURLOPT_MAXREDIRS => 20,
        // Callers check for "HTTP/1.1 200" in the status line
        CURLOPT_HTTP_VERSION => CURL_HTTP_VERSION_1_1,
        CURLOPT_HTTPHEADER => $request_headers,
        CURLOPT_CONNECTTIMEOUT => 30,
        CURLOPT_TIMEOUT => (int) ($http['timeout'] ?? (ini_get('default_socket_timeout') ?: 60)),
        CURLOPT_HEADERFUNCTION => function ($ch, $line) use (&$response_headers) {
            $trimmed = trim($line);
            // A new status line means a new response (after a redirect): only keep the final one
            if (preg_match('#^HTTP/\S+\s+\d{3}#', $trimmed)) {
                $response_headers = [];
            }
            if ($trimmed !== '') {
                $response_headers[] = $trimmed;
            }
            return strlen($line);
        },
    ]);
    if ($method === 'POST') {
        curl_setopt($ch, CURLOPT_POST, true);
    } elseif ($method !== 'GET') {
        curl_setopt($ch, CURLOPT_CUSTOMREQUEST, $method);
    }
    if (isset($http['content']) && $method !== 'GET') {
        curl_setopt($ch, CURLOPT_POSTFIELDS, $http['content']);
    }
    if (isset($http['user_agent'])) {
        curl_setopt($ch, CURLOPT_USERAGENT, $http['user_agent']);
    }

    $body = curl_exec($ch);
    $status = curl_getinfo($ch, CURLINFO_RESPONSE_CODE);
    $error = curl_error($ch);
    fetch_url_headers($response_headers);

    if ($body === false) {
        trigger_error("fetch_url($url): $error", E_USER_WARNING);
        return false;
    }
    if ($status >= 400 && empty($http['ignore_errors'])) {
        trigger_error(
            "fetch_url($url): HTTP request failed! " . ($response_headers[0] ?? "HTTP status $status"),
            E_USER_WARNING,
        );
        return false;
    }
    return $body;
}
