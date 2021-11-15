
import requests
import time
import re
from Bio.Blast.NCBIWWW import qblast


class BlastApi:
    @staticmethod
    def send_query(program: str, database: str, seq: str) -> bool:
        """
        Methode which connects to blast.api and saves results from queried seq. to a file.
        :param program: blastn, blastp, blastx, tblastn, tblastx
        :param database: nr, cdd, nt, ...
        :param seq: sequence of nucleotids ACGT
        :return: True/False if the query retrieved some hits or not
        """
        try:
            program = "blastn"
            database = "nt"
            seq = "AGGATATTGTATTAGACCTGCAACCTCCAGACCCTGTAGGGTTACATTGCTATGAGCAATTAGTAGACAGCGCAGA"
            sequence_data = open("blast_example.fasta").read()

            url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params = {"CMD" : "Put", "PROGRAM" : "blastn", "DATABASE" : "nt", "QUERY" : seq}
            headers = {'Content-type': 'application/x-www-form-urlencoded'}
            query_request = requests.post(url=url_base, headers=headers, params=params)

            print(query_request.status_code)

            # parse out the estimated time to completion
            index_rtoe = query_request.content.index(str.encode("RTOE ="))
            rtoe = int(query_request.content[index_rtoe:].decode().split("\n")[0].split("=")[1].lstrip())

            # parse out the request id
            index_rid = query_request.content.index(str.encode("RID ="))
            rid = query_request.content[index_rid:].decode().split("\n")[0].split("=")[1].lstrip()

            start = time.time()

            # wait for search to complete
            time.sleep(rtoe)

            # poll for results
            while True:
                time.sleep(5)

                url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
                params.clear()
                params = {"CMD" : "Get", "FORMAT_OBJECT" : "SearchInfo", "RID" : rid}
                result_ready_request = requests.get(url=url_base, params=params)
                print(result_ready_request.status_code)

                index_status = result_ready_request.content.index(str.encode("Status="))
                status = result_ready_request.content[index_status:].decode().split("\n")[0].split("=")[1].lstrip()

                if status == "WAITING":
                    print("Searching...")
                    continue
                elif status == "FAILED":
                    end = time.time()
                    print(f"Time of search: {end - start}")
                    raise SystemError(f"Search {rid} failed; please report to blast-help\\@ncbi.nlm.nih.gov.")
                elif status == "UNKNOWN":
                    end = time.time()
                    print(f"Time of search: {end - start}")
                    raise SystemError(f"Search {rid} expired.")
                elif status == "READY":
                    end = time.time()
                    print(f"Time of search: {end-start}")

                    index_hits = result_ready_request.content.index(str.encode("ThereAreHits="))
                    there_are_hits = result_ready_request.content[index_hits:].decode().split("\n")[0].split("=")[1].lstrip()
                    if there_are_hits == "yes":
                        print("Search is complete, retrieving results...")
                        break
                    else:
                        print("No hits found.")
                        return False
                else:
                    end = time.time()
                    print(f"Time of search: {end - start}")
                    raise SystemError("Something unexpected happened during search. Please try it later.")

            # retrieve and display results
            url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params.clear()
            params = {"CMD": "Get", "FORMAT_TYPE": "Text", "RID": rid}
            result_request = requests.get(url=url_base, params=params)
            print(result_request.status_code)

            with open('results.txt', 'w') as save_file:
                save_file.write(result_request.content.decode())

            return True

        except Exception as e:
            print(f"Oops, something bad happened : {repr(e)}")
            return False









