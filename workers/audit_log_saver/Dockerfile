# Use official Python image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install Poetry
RUN pip install poetry==1.8.3

# Copy message_queue into the container.
COPY --from=webapp /sources/services/message_queue /app/Webapp/sources/services/message_queue

# Copy the 'audit_log_saver' directory into the container.
COPY . /app/workers/audit_log_saver

ENV PYTHONPATH="/app/Webapp"
WORKDIR /app/workers/audit_log_saver

# Install dependencies via Poetry
RUN poetry install --no-dev

# Command to run the worker script
CMD ["poetry", "run", "python", "audit_log_saver.py"]
