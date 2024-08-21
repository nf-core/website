-- FIXME This file isn't creating a "users_data" table
-- I guess you can't write sql queries for local csvs?
SELECT
    to_timestamp(unix_time) AS timestamp,
    total_users,
    active_users,
    inactive_users
FROM slack.users
