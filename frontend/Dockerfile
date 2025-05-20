# ── build stage ─────────────────────────────────────────────
FROM node:20-alpine AS build
WORKDIR /src
COPY frontend/package*.json ./
RUN npm ci --legacy-peer-deps
COPY frontend .
RUN npm run build      # produces static files in /src/public

# ── runtime image ───────────────────────────────────────────
FROM nginx:1.25-alpine
# Replace default Nginx site with our app
COPY --from=build /src/public /usr/share/nginx/html
COPY frontend/nginx.conf /etc/nginx/conf.d/default.conf
# Optional: tailor cache rules / gzip, etc.
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"] 