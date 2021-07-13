FROM node:15.2.0-alpine3.12

WORKDIR .

COPY package.json .

RUN mv /usr/local/lib/node_modules /usr/local/lib/node_modules.tmp \
    && mv /usr/local/lib/node_modules.tmp /usr/local/lib/node_modules && \
    npm install -g npm && npm install

COPY node_modules .
