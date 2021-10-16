import React, { Component } from 'react';
import { BrowserRouter as Router, Route } from "react-router-dom";
import Board from './board';
import Instructions from './instructions';
import LeaderBoard from './leaderboard';
import WelcomeScreen from './welcomeScreen';
import Header from './header';
import firebase from '../firebase.js';
import queryString from 'query-string';
import './app.css';
import '../tachyons.css';

class App extends Component {
  constructor(props) {
    super(props);

    this.state = {
      gameId: queryString.parse(props.location.search).game,
      noSuchGame: false,
      signedIn: false,
    };

    firebase.database().ref('games/' + this.state.gameId).once('value').then((game) => {
      this.setState({
        noSuchGame: !game.exists(),
        lexicon: game.child('lexicon').val(),
        size: game.child('size').val()
      });
      let instructions = game.child('instructions').val();
      if (instructions) {
        this.setState({
          instructions: instructions.replace(/\\n/g, '\n')
        });
      }
    });

    firebase.database().ref('games/' + this.state.gameId + '/leaderboard').orderByChild('duration').on('value', (leaderboard) => {
      const values = leaderboard.val();
      const leaders = 
        values ?
          Object.values(leaderboard.val()).sort((a, b) => {
            return a.duration - b.duration;
          }) : [];
      this.setState({leaders: leaders});
    });
  }

  renderGame() {
    return (
      <div>
        <Board id='abc' size={this.state.size} values={this.state.lexicon} db={firebase} gameId={this.state.gameId} />
        <LeaderBoard leaders={this.state.leaders} size={this.state.size} />
        <Instructions src={this.state.instructions}/>
      </div>
    );
  }

  renderLoadingScreen() {
    return (
      <div>
        <Header/>
        <main>
          <div className="f3 pa2" aria-live='polite'>Loading...</div>
        </main>
      </div>
    );
  }

  renderWelcomeScreen() {
    return (
      <div>
        <Header/>
        <main>
          <WelcomeScreen firebase={firebase} />
        </main>
        <Instructions/>
      </div>
    );
  }

  render() {
    if (this.state.gameId) {
      if (!this.state.noSuchGame) {
        if (this.state.lexicon && this.state.size) {
          return this.renderGame();
        } else {
          return this.renderLoadingScreen();
        }
      } 
    }
    return this.renderWelcomeScreen();
  }
}

const AppRouter = () => (
  <Router>
    <Route path="/" exact component={App} />
  </Router>
);

export default AppRouter;
