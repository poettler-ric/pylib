#!/usr/bin/env python3

from statistics import mean

from berserk import Client


def get_opponents(client, player, max_games=20):
    games = client.games.export_by_player(player, max=max_games)
    opponents = []
    score = 0
    for game in games:
        white_player = game["players"]["white"]["user"]["name"]
        black_player = game["players"]["black"]["user"]["name"]
        if game["status"] == "draw":
            score += 0.5
        elif game["winner"] == "white" and white_player == player:
            score += 1
        opponents.append(white_player if white_player != player else black_player)

    return opponents, score


def get_game_count(client, player, controls=["blitz", "rapid"]):
    count = []
    performance = client.users.get_public_data(player)
    for c in controls:
        count.append(performance["perfs"][c]["games"])
    return count


def main():
    __USERNAME = "Iwirada"
    __GAMES_PER_DAY = 4
    __WEEKS = 4

    client = Client()

    game_count = sum(get_game_count(client, __USERNAME))
    print(f"{__USERNAME} game count: {game_count}")

    opponents, _ = get_opponents(client, __USERNAME, __GAMES_PER_DAY * __WEEKS * 7)
    game_counts = []
    for opponent in opponents:
        game_counts.append(sum(get_game_count(client, opponent)))
    print(f"mean opponent game count: {mean(game_counts):.1f}")


if __name__ == "__main__":
    main()
