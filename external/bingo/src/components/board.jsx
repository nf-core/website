import React, { Component } from 'react';
import PropTypes from 'prop-types';
import moment from 'moment';
import momentDurationFormatSetup from 'moment-duration-format';
import SimpleStorage from 'react-simple-storage';
import Header from './header';

momentDurationFormatSetup(moment);

function copyObject(obj) {
  return Object.assign({}, obj);
}

/* https://stackoverflow.com/a/12646864 */
function shuffleArray(array) {
  for (let i = array.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [array[i], array[j]] = [array[j], array[i]];
  }

  return array;
}

class Board extends Component {

  constructor(props) {
    super(props);

    const size = props.size % 2 ? props.size : props.size - 1;
    const cellCount = size * size;
    const midpoint = (size * size - 1)/ 2;

    // Ensure we have enough values to fill this size board.
    // If not, keep adding duplicate values until we do.
    let values = props.values.slice();
    let i = 0;
    while (values.length < cellCount) {
      values.push(props.values[i]);
      i++;
      if (i > props.values.length - 1) i = 0;
    }

    this.state = {
      activeCell: 0,
      activeRow: 0,
      activeCol: 0,
      endTime: 0,
      grid: this.generateRandomGrid(values, size),
      midpoint: midpoint,
      selection: {[midpoint]: true},
      size: size,
      startTime: Date.now(),
      values: values,
    };

    this.handleKeyDown = this.handleKeyDown.bind(this);
    this.refreshBoard = this.refreshBoard.bind(this);
    this.updateLeaderBoard = this.updateLeaderBoard.bind(this);
  }

  /**
   * Randomize supplied values and return 
   * a grid with dimensions size * size
   */
  generateRandomGrid(values, size) {
    const randomizedValues = shuffleArray(values);

    let grid = [];
    for (let row = 0; row < size; row++) {
      grid[row] = [];
      for (let col = 0; col < size; col++) {
        let id = col + (row * size);
        grid[row][col] = {
          value: randomizedValues[id],
          id: id
        }
      }
    }

    return grid;
  }

  updateLeaderBoard() {
    const name = document.getElementById('name').value;
    if (name && name !== '') {
      const data = {
        name: name,
        timestamp: this.state.endTime,
        duration: this.state.endTime - this.state.startTime
      };
      const key = this.props.db.database().ref('games/' + this.props.gameId).child('leaderboard').push(data).key;
      this.setState({leaderboardSubmitted: true});
      return key;
    }
  }

  /**
   * Randomize cell values, reset timer, and clear selection.
   */
  refreshBoard() {
    this.setState({
      activeCell: 0,
      activeRow: 0,
      activeCol: 0,
      bingo: false,
      grid: this.generateRandomGrid(this.state.values, this.state.size),
      leaderboardSubmitted: false,
      selection: {[this.state.midpoint]: true},
      startTime: Date.now(),
      endTime: 0
    });
  }

  componentDidUpdate(prevProps, prevState) {
    // not a new board
    if (prevState.startTime === this.state.startTime) {
      // focus active cell
      if (prevState.activeCell !== this.state.activeCell) {
        document.getElementById(this.props.id + '-cell-' + this.state.activeCell).focus();
      }

      // if selection has changed in some way, check for bingo
      if (prevState.selection !== this.state.selection) {
        console.clear();
        console.log('Checking for bingo...');
        if (
          this.checkRow(this.state.activeRow) ||
          this.checkCol(this.state.activeCol) ||
          this.checkDiagonalA(this.state.activeRow, this.state.activeCol) ||
          this.checkDiagonalB(this.state.activeRow, this.state.activeCol)
        ) {
          if (!this.state.bingo) {
            console.log('It\'s a bingo!');
            this.setState({
              bingo: true,
              endTime: Date.now(),
            });
          }
        } else {
          console.log('Alas, not a bingo.')
          this.setState({bingo: false});
        }
      }
    }
  }

  checkIndices(indices) {
    for (let i = 0; i < indices.length; i++) {
      let index = indices[i];
      if (!this.state.selection[index]) {
        return false;
      }
    }
    return true;
  }

  checkRow(row) {
    const size = this.state.size;
    const rowStart = row * size;
    for (let i = rowStart; i < rowStart + size; i++) {
      if (!this.state.selection[i]) {
        console.log('- Bingo in row ' + row + '? Fail at cell ' + i);
        return false;
      }
    }
    console.log('- Bingo in row ' + row + '? Success!');
    return true;
  }

  checkCol(col) {
    const size = this.state.size;
    for (let j = col; j < size * size; j+= size) {
      if (!this.state.selection[j]) {
        console.log('- Bingo in col ' + col + '? Fail at cell ' + j);
        return false;
      }
    }
    console.log('- Bingo in col ' + col + '? Success!');
    return true;
  }

  /* Upper left to lower right */
  checkDiagonalA(row, col) {
    const size = this.state.size;
    if (row === col) {
      for (let i = 0; i < size; i++) {
        if (!this.state.selection[size * i + i]) {
          console.log('- Bingo in diagonal A? Fail at cell ' + (size * i + i));
          return false;
        }
      }
      console.log('- Bingo in diagonal A? Success!');
      return true;
    }
  }

  /* Upper right to lower left */
  checkDiagonalB(row, col) {
    const size = this.state.size;
    if (row === (size - col - 1)) {
      for (let i = 0; i < size; i++) {
        if (!this.state.selection[size * i + size - i - 1]) {
          console.log('- Bingo in diagonal B? Fail at cell ' + (size * i + size - i - 1));
          return false;
        }
      }
      console.log('- Bingo in diagonal B? Success!')
      return true;
    }
  }

  handleKeyDown(event, row, col) {
    switch (event.key) {
      case 'Down':
      case 'ArrowDown':
        if (row < this.state.size - 1) this.setActiveCell(row + 1, col);
        event.preventDefault();
        break;
      case 'Up':
      case 'ArrowUp':
        if (row > 0) this.setActiveCell(row - 1, col);
        event.preventDefault();
        break;
      case 'Left':
      case 'ArrowLeft':
        if (col > 0) this.setActiveCell(row, col - 1);
        event.preventDefault();
        break;
      case 'Right':
      case 'ArrowRight':
        if (col < this.state.size - 1) this.setActiveCell(row, col + 1);
        event.preventDefault();
        break;
      default:
        break;
    }
  }

  setActiveCell(row, col) {
    this.setState({activeCell: this.state.grid[row][col].id});
  }

  renderMidpointCell(cellId, row, col) {
    return (
      <td role='gridcell' key={cellId}>
        <div className='cell-contents'>
          <button
            aria-disabled={true}
            aria-pressed={true}
            className='cell-toggle'
            id={this.props.id + '-cell-' + cellId}
            onClick={() => {this.setState({activeCell : cellId});}}
            onKeyDown={(event) => {this.handleKeyDown(event, row, col);}}
            tabIndex={cellId === this.state.activeCell ? '0' : '-1'}
          >
            <svg role="img" aria-labelledby="star-title" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
              <title id="star-title">Star (free tile)</title>
              <path d="M12.6 1.4l2.2 7c.1.2.3.4.6.4h6.9c.7 0 1 .9.5 1.3l-5.7 4.2c-.2.1-.3.5-.2.7l2.7 7.2c.2.6-.5 1.2-1.1.7l-6-4.5c-.3-.2-.6-.2-.9 0l-6.1 4.5c-.5.5-1.3-.1-1-.7L7.1 15c.1-.2 0-.6-.3-.7l-5.6-4.2c-.6-.4-.2-1.3.4-1.3h6.9c.4 0 .6-.1.7-.4l2.2-7c.1-.7 1.1-.6 1.2 0z"></path>
            </svg>
          </button>
        </div>
      </td>
    );
  }

  renderCell(cell, row, col) {
    const isMidpoint = cell.id === this.state.midpoint;
    const selected = this.state.selection[cell.id] || isMidpoint ? true : false;

    if (isMidpoint) { return this.renderMidpointCell(cell.id, row, col); }

    return (
      <td role='gridcell' key={cell.id}>
        <div className='cell-contents'>
          <button
            aria-pressed={selected}
            className='cell-toggle'
            id={this.props.id + '-cell-' + cell.id}
            onClick={() => {
              let selection = copyObject(this.state.selection);
              selection[cell.id] = !selected;

              this.setState({
                selection: selection,
                activeCell: cell.id,
                activeRow: row,
                activeCol: col
              });
            }}
            onKeyDown={(event) => {this.handleKeyDown(event, row, col);}}
            tabIndex={cell.id === this.state.activeCell ? '0' : '-1'}
          >
            {/*cell.id*/}
            {cell.value} 
          </button>
        </div>
      </td>
    );
  }

  renderRow(row, y) {
    return (
      <tr role='row' key={y}>
        {row.map((cell, x) => { return this.renderCell(cell, y, x); })}
      </tr>
    );
  }

  renderLeaderboardPrompt() {
    if (this.state.bingo) {
      if (this.state.leaderboardSubmitted) {
        return (
          <p className='lh-copy mb0'>
            You're on the leaderboard!{' '}
            Keep playing on this bingo board or{' '}
            <button onClick={this.refreshBoard}>generate a new one</button>.
          </p>
        );
      } else {
        return (
          <div className='pt3'>
            <label htmlFor='name'>Enter a name to display on the leaderboard:</label>
            <div className='pv2'>
              <input
                style={{'backgroundColor': '#f6f7fa'}}
                className='input-reset pa3 mb2 ba bw1 b--black'
                id='name'
                placeholder='Your Name' />
              <button
                className='tc fw8 blue-button white pa3 ml2 ba bw1 b--black'
                onClick={this.updateLeaderBoard}
              >
                Add me!
              </button>
            </div>
          </div>
        );
      }
    }
    return null;
  }

  renderSuccess() {
    if (this.state.bingo) {
      return (
        <div className='success maxw-95 pa3 mv3'>
          <div className='flex flex-wrap items-center justify-between'>
            <div className='w-50-l w-100 tc tl-l' role='alert' aria-live='assertive'>
              <span className='f2 fw8'>You got bingo! <span role="img" aria-label="Hurray!">🎉</span></span>
              <div className='f3 pt2'>
                Total time: {moment.duration(this.state.endTime - this.state.startTime).format('h [hr], m [min], s [sec]')}
              </div>
            </div>
            <div className='w-50-l w-100 tc tr-l'>
              {this.renderLeaderboardPrompt()}
            </div>
          </div>
        </div>
      )
    }
    return null;
  }

  render() {
    return (
      <div>
        <Header gotBingo={this.state.bingo}>
          <button
            className='tc fw8 bg-white black pa3 ba bw1 b--black mb2'
            onClick={this.refreshBoard}
          >
            New Board
          </button>
        </Header>
        <main>
          <p style={{fontSize: "12px"}}>🛈 See leaderboard and Instructions below the board.</p>
          <table role='grid' className='maxw-95'>
            <tbody role='rowgroup'>
              {this.state.grid.map((row, y) => { return (this.renderRow(row, y))})}
            </tbody>
          </table>
          {this.renderSuccess()}
        </main>
        { /* Stores current board state in local storage so
             game is preserved even when refreshed */ }
        <SimpleStorage
          parent={this}
          prefix={`bingo-${this.props.gameId}`}
          blacklist={['activeCell', 'activeRow', 'activeCol', 'values']}
        />
      </div>
    );
  }
}

Board.propTypes = {
  size: PropTypes.number,
  values: PropTypes.array
}

Board.defaultProps = {
  size: 5,
  values: 'abcdefghijklmnopqrstuv'.split('')
}

export default Board;
