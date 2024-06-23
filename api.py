import time
import logging
import colorlog
from typing import Dict, Any

from fastapi import FastAPI, Request, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from utils.processor import Field

handler = colorlog.StreamHandler()
handler.setFormatter(colorlog.ColoredFormatter(
    "%(log_color)s%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    log_colors={
        'DEBUG': 'white',
        'INFO': 'green',
        'WARNING': 'yellow',
        'ERROR': 'red',
        'CRITICAL': 'red,bg_white',
    }))

# Configure logging to write logs to a file
file_handler = logging.FileHandler("app.log")
file_handler.setFormatter(colorlog.ColoredFormatter(
    "%(log_color)s%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    log_colors={
        'DEBUG': 'white',
        'INFO': 'green',
        'WARNING': 'yellow',
        'ERROR': 'red',
        'CRITICAL': 'red,bg_white',
    }))


# Create a logger and add the handlers
logger = logging.getLogger()
logger.addHandler(handler)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO) 

"""
uvicorn api:app --reload --port 8000
"""

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

logger.info("Starting the API")

field = Field("/home/cesaire/Projects/terraventis/farmland/6.geojson")

@app.get("/indices", response_model=Dict[str, Any])
async def get_indices(start_date: str, end_date:str):
    try:
        aggregated_data = field.computeIndexes(start_date, end_date)
        return {"response": aggregated_data}
    except:
        raise HTTPException(status_code=404, detail="Date not found")


# Middleware for logging and request timing
@app.middleware("http")
async def log_requests(request: Request, call_next):
    '''
    Middleware for logging and request timing
    '''
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    logger.info(f"{request.method} {request.url.path} - Process Time: {process_time:.2f}")
    return response

# Endpoint to display logs
@app.get("/logs", response_class=HTMLResponse)
async def display_logs():
    '''
    Endpoint to display logs
    '''
    logger.info("Displaying logs")
    # Open and read the log file
    with open("app.log", "r") as log_file:
        logs = log_file.readlines()
        logs_html = "<br>".join(logs)
        return f"<h1>Logs</h1><pre>{logs_html}</pre>"
