# Experimental - not being used right now

import asyncio
import websockets
import threading

PORT = 8034


def ws_server(cmd_pointer):
    # workspace = cmd_pointer.settings["workspace"]
    # data = {"foo": 123, "log": []}

    # print(f"Starting websocket server...\nws://localhost:{8025}")
    # print("% workspace:", workspace)

    async def hello(websocket, path):
        message = await websocket.recv()
        # print(f"Received message: {message}")
        # data["log"].append(message)
        # print(data)

        response = "Hello, World!"
        await websocket.send(response)
        # print(f"Sent message: {response}")

    def start_server():
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        start_server = websockets.serve(hello, "127.0.0.1", PORT)

        loop.run_until_complete(start_server)
        loop.run_forever()
        print("Started")

    # Start the WebSocket server in a new thread
    threading.Thread(target=start_server).start()


# def send(message):
#     print(f"Sending message: {message}")
#     websocket.send(response)
#     return message


async def send(message):
    uri = f"ws://127.0.0.1:{PORT}"
    async with websockets.connect(uri) as websocket:
        await websocket.send(message)
        # print(f"Message sent: {message}")
