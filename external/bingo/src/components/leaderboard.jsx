import React, { Component } from 'react';
import moment from 'moment';
import momentDurationFormatSetup from 'moment-duration-format';

momentDurationFormatSetup(moment);

class LeaderBoard extends Component {

  renderEmptyState() {
    if (!this.props.leaders || !this.props.leaders.length) {
      return (
        <div>
          <p>There's no one on the leaderboard yet!</p>
          <p>Get {this.props.size} phrases in a row to add yourself to the list with: &lt;name surname&gt; (@&lt;github handle&gt;).</p>
        </div>
      );
    }
    return null;
  }

  render() {
    return (
      <aside className='maxw-95' aria-label="Leaderboard">
        <h2 className="bb-3 pv2">Leaderboard</h2>
        <div role='log' aria-live='polite' aria-atomic='true'>
          <ol className="f4 list pa0">
            {this.props.leaders.map((leader, i) => {
              return (
                <li key={i} className='ph1 pb3'>
                  <span aria-hidden="true">🏆 </span> 
                  <strong>
                    {leader.name}
                    {' '}&middot;{' '}
                    {moment.duration(leader.duration).format('h [hr], m [min], s [sec]')}
                    {' '}&middot;{' '}
                  </strong>
                  {moment(leader.timestamp).format('L LT')}
                </li>
              );
            })}
          </ol>
          {this.renderEmptyState()}
        </div>
      </aside>
    );
  }

}

export default LeaderBoard;