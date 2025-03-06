FROM node:${NODE_VERSION}

WORKDIR /var/www/html/app
USER node

COPY --chown=node:node package*.json ./

RUN npm install

COPY --chown=node:node . .

EXPOSE 5173