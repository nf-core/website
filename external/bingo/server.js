const express = require('express');
const http = require('http');
const socketIO = require('socket.io');

const port = 4001;

const app = express();

const server = http.createServer(app);

const io = socketIO(server);

io.on('connection', socket => {
  console.log('User connected')
  
  socket.on('disconnect', () => {
    console.log('user disconnected')
  })
})

server.listen(port, () => console.log(`Listening on port ${port}`))