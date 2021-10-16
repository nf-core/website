# Bingo

An accessible bingo web app using React and Firebase. Modified from the 'backlog-bingo' app by [Cordelia Dillon](https://github.com/cordeliadillon) for the March 2021 nf-core/hackathon.

> All documentation below is from the original author

Original:
[See it in action](https://backlog-bingo.com/), [read about its accessibility](https://www.24a11y.com/2019/building-an-accessible-bingo-web-app/), or build your own version using the instructions below:

## Getting started

### Database setup

You could use any underlying database for this, but [Firebase](https://firebase.google.com/) is easy and has nice realtime update support.

1. Create a new Firebase project at https://firebase.google.com/.
2. In the Firebase console for your project, create a new Realtime Database.
3. Specify location in US and start in locked mode.
4. In the Data tab of your Realtime Database, click Import JSON and upload the [sample schema](sample_schema.json).
5. Copy the rules in `database_rules.json` into the Realtime Database > Rules tab.

> N.B. Make sure Default GCP resource location is US (Europe will not work, as this doesn't give the `firebase.io` URL suffix)

### Local setup

1. Clone this repo.
2. Run `npm install`.
3. Run `firebase init`. Follow the setup steps, making sure to:
    - choose Firebase CLI features **Database** & **Hosting**
    - choose the Firebase project you just created
    - not overwrite `database.rules.json`
    - choose `build` as your public directory
    - configure as single-page app
4. Create an `.env` file with the following data from your Firebase project's settings:
    ```
    REACT_APP_FIREBASE_API_KEY={your_data} ## Need to activate by loading Build > Authentication tab (and turn on anonymous at minimum). Find API key under project settings.
    REACT_APP_FIREBASE_AUTH_DOMAIN={your_data} ## Not Needed
    REACT_APP_FIREBASE_DATABASE_URL={your_data} #  Build > Realtime Database  > As displayed on the link above the JSON thing on
    REACT_APP_FIREBASE_PROJECT_ID={your_data} ## Project Settings > General > Project ID
    REACT_APP_FIREBASE_STORAGE_BUCKET={your_data} # ## Project Settings > General > Default GCP resource location (Must be US as URL needs to end in firebase.io)
    REACT_APP_FIREBASE_MESSAGEING_SENDER_ID={your_data} ## Project Settings > General > Project Number
    ```

### Check that it's all working

If you've followed the steps above, you should be able to run the app locally.

1. Call `npm start`.
2. App should automatically load showing the "Ready to join a game?" screen.
3. If you uploaded the sample schema, type "a11y" into the Board Name field and click "Play!"

If the app hangs at this point, go into your Firebase Console and make sure your database's rules are the same as those in [database.rules.json](database.rules.json).


## Creating new games

Right now this is a manual process as I've yet to build the authenticated UI for doing this within the bingo app itself. Whoops!

The easiest way to add a new game board is to create a JSON file containing a lexicon of phrases, a numerical size (how many squares across the bingo grid will be), and optionally a set of instructions for your users. Like so:

```
{
  "instructions" : "### Write anything here using basic **markdown**...",
  "lexicon" : [
    "Put all the phrases you want to use here",
    "You can even include emojis!",
    "🍕"
   ],
  "size" : 5
}

```

In your Firebase Console, create a new entry in the "games" array and import your JSON into that object. Whatever name you choose for that entry will be the "Board Name" players will use to access your game.


## Deploying

Should be as simple as:

```
> npm run build
> firebase deploy
```

> If you get an error such as `Error: HTTP Error: 403, Cloud Firestore API has not been used in project`, you need to enable the Firestore API. Use the printed link, but check you open with the right Google account.

> If you get an error such as `Error: HTTP Error: 400, <...something about doesn't work with datastore anymore>`, go to https://console.cloud.google.com/firestore/ and activate `Native` mode, and try deploying again

> You are looking for 15 files in `build` in total being deployed

## Disclaimers

This is a prototype for which I have yet to build any tests. Use at your own risk.

## Credits

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).
