#!/usr/bin/env python3
"""
---
requirements:
  geopy: ">=2.4.1"
---
Script to generate location markdown files from CSV responses for the nf-core hackathon.
"""

import csv
import re
import os
from geopy.geocoders import Nominatim
from geopy.exc import GeocoderTimedOut, GeocoderUnavailable
import time


def get_coordinates(address, org_name):
    geolocator = Nominatim(user_agent="nf-core-hackathon-sites")
    try:
        # Try with full address first
        location = geolocator.geocode(address)
        if location:
            return [location.latitude, location.longitude]

        # If that fails, try with organization name
        location = geolocator.geocode(org_name)
        if location:
            return [location.latitude, location.longitude]

        # If both fail, return None
        return None
    except (GeocoderTimedOut, GeocoderUnavailable):
        # Handle timeout gracefully
        print(f"Warning: Geocoding failed for {org_name}")
        return None
    finally:
        # Be nice to the geocoding service
        time.sleep(1)


def slugify(text):
    # Custom mappings for specific organizations
    mappings = {
        "brain-institute-federal-university-of-rio-grande-do-norte": "university-of-rio-grande-do-norte",
        "dept-biology-and-biotechnology-university-of-pavia": "university-of-pavia",
        "dept-chemistry-and-biology-a-zambelli-university-of-salerno": "university-of-salerno",
        "the-german-cancer-research-center-and-ghga": "dkfz-and-ghga"
    }

    # First do the normal slugification
    text = text.lower().strip()
    # Remove special characters and replace spaces/underscores with hyphens
    text = re.sub(r'[^\w\s-]', '', text)
    text = re.sub(r'[-\s]+', '-', text)
    text = text.replace('ü', 'u').replace('ä', 'a').replace('ö', 'o')

    # Then check if we need to override with a custom mapping
    return mappings.get(text, text)


def parse_timezone(tz_string):
    # Extract GMT offset from timezone string
    match = re.search(r'GMT([+-])(\d+)', tz_string)
    if match:
        sign = match.group(1)
        hours = int(match.group(2))
        return f"{sign}{hours:02d}:00"
    return "+00:00"  # Default to UTC if no match


def create_markdown_file(row, output_dir):
    # Extract organization name and create slug
    org_name = row["What's the name of the organization where the hackathon will take place?"]
    location_slug = slugify(org_name)

    # Parse timezone
    tz_offset = parse_timezone(row["What time zone is your site in?"])

    # Get coordinates
    venue = row["What's the venue that's been booked? "]
    coordinates = get_coordinates(venue, org_name)
    if coordinates is None:
        coordinates = [0, 0]
        print(f"Warning: Using default coordinates for {org_name}")

    # Create frontmatter content
    frontmatter = {
        "title": f"Hackathon - March 2025 ({org_name})",
        "subtitle": f"Local node of the nf-core hackathon at {org_name}",
        "shortTitle": org_name,
        "type": "hackathon",
        "startDate": "2025-03-17",  # Hardcoded dates for March 2025
        "startTime": f"09:00{tz_offset}",
        "endDate": "2025-03-19",
        "endTime": f"17:00{tz_offset}",
        "locations": [{
            "name": org_name,
            "address": venue,
            "links": [row["What's {{field:82459451-2436-4a71-8c6b-d2e64e3f2f6c}}'s website?"]],
            "geoCoordinates": coordinates
        }]
    }

    # Create markdown content
    content = "---\n"
    for key, value in frontmatter.items():
        if key == "locations":
            content += "locations:\n"
            for location in value:
                content += "  - name: " + str(location["name"]) + "\n"
                content += "    address: |\n      " + str(location["address"]).replace("\n", "\n      ") + "\n"
                content += "    links:\n"
                for link in location["links"]:
                    content += f"      - {link}\n"
                content += "    geoCoordinates: " + str(location["geoCoordinates"]) + "\n"
        else:
            content += f"{key}: '{str(value)}'\n"
    content += "---\n\n"  # Add extra newline after frontmatter

    # Add body content
    org_name = row["What's the name of the organization where the hackathon will take place?"]
    contact_name = f"{row['First name']} {row['Last name']}"
    contact_email = row['Email']
    external_allowed = row["Will your site be open to attendees from outside of {{field:82459451-2436-4a71-8c6b-d2e64e3f2f6c}}?"]
    venue_booked = row["Have you already booked a space for hackathon?"]

    # Add attendance scope text
    if external_allowed == "1":
        content += f"Attendees from outside {org_name} **can** attend.\n\n"
    else:
        content += f"Attendees from outside {org_name} **cannot** attend.\n\n"

    # Add contact information
    content += f"Main contact: {contact_name} ([{contact_email}](mailto:{contact_email}))\n\n"

    # Add venue booking status
    if venue_booked != "Yes":
        content += "Venue not yet booked, please check back later for more information.\n\n"

    # Write to file
    filename = f"{location_slug}.md"
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

    return filename


def main():
    # Get the directory of the CSV file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(current_dir, "responses.csv")

    # Read CSV and create markdown files
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            filename = create_markdown_file(row, current_dir)
            print(f"Created {filename}")


if __name__ == "__main__":
    main()
