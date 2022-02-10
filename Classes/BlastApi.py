import threading

import requests
import time
import re
from Bio.Blast.NCBIWWW import qblast

from Classes.BlastResult import BlastResult
from Classes.SamRecord import SamRecord


blast_results = []

states = ["NO_HIT", "UNKNOWN", "FAILED", "UNEXPECTED", "SUCCESS"]

class BlastApi:
    @staticmethod
    def send_multiple_queries(program: str, database:str, sam_records: list[SamRecord]):
        blast_results.clear()
        threads = []

        start = time.time()
        for no, sam_record in enumerate(sam_records):
            t = threading.Thread(target=BlastApi.send_query, args=(program, database, sam_record, no,))
            threads.append(t)
            t.start()
        for t in threads:
            t.join()

        end = time.time()
        print(f"Time of Blast Api work: {end - start}")
        return blast_results
    @staticmethod
    def send_query(program: str, database: str, sam_record: SamRecord, no: int) -> bool:
        """
        Methode which connects to blast.api and saves results from queried seq. to a file.
        :param no: number of record.
        :param program: blastn, blastp, blastx, tblastn, tblastx
        :param database: nr, cdd, nt, ...
        :param sam_record: Record that contain sequence of nucleotids ACGT
        :return: True/False if the query retrieved some hits or not
        """
        try:
            seq = sam_record.SEQ
            url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params = {"CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": seq}
            headers = {'Content-type': 'application/x-www-form-urlencoded'}
            query_request = requests.post(url=url_base, headers=headers, params=params)

            print(f"[Record no. {no}.]:{query_request.status_code}")

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
                time.sleep(60)

                url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
                params.clear()
                params = {"CMD" : "Get", "FORMAT_OBJECT" : "SearchInfo", "RID" : rid}
                result_ready_request = requests.get(url=url_base, params=params)


                index_status = result_ready_request.content.index(str.encode("Status="))
                status = result_ready_request.content[index_status:].decode().split("\n")[0].split("=")[1].lstrip()
                print(f"[Record no. {no}.]:{result_ready_request.status_code}, status: {status}")

                if status == "WAITING":
                    print(f"[Record no. {no}.]: Searching...")
                    continue
                elif status == "FAILED":
                    end = time.time()
                    print(f"[Record no. {no}.]: FAILED. Time of search: {end - start}")
                    blast_results.append(BlastResult(sam_record, "FAILED", rid, no, None))
                    raise SystemError(f"[Record no. {no}.]: Search {rid} failed; please report to blast-help\\@ncbi.nlm.nih.gov.")
                elif status == "UNKNOWN":
                    end = time.time()
                    print(f"[Record no. {no}.]: UNKNOWN. Time of search: {end - start}")
                    blast_results.append(BlastResult(sam_record, "UNKNOWN", rid, no, None))
                    raise SystemError(f"[Record no. {no}.]: Search {rid} expired.")
                elif status == "READY":
                    end = time.time()
                    print(f"[Record no. {no}.]: READY. Time of search: {end-start}")

                    index_hits = result_ready_request.content.index(str.encode("ThereAreHits="))
                    there_are_hits = result_ready_request.content[index_hits:].decode().split("\n")[0].split("=")[1].lstrip()
                    if there_are_hits == "yes":
                        print(f"[Record no. {no}.]: Search is complete, retrieving results...")
                        break
                    else:
                        print(f"[Record no. {no}.]: No hits found.")
                        blast_results.append(BlastResult(sam_record, "NO_HITS", rid, no, None))
                        return False
                else:
                    end = time.time()
                    print(f"[Record no. {no}.]: UNEXPECTED. Time of search: {end - start}")
                    blast_results.append(BlastResult(sam_record, "UNEXPECTED", rid, no, None))
                    raise SystemError(f"[Record no. {no}.]: Something unexpected happened during search. Please try it later.")

            # retrieve and display results
            url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
            params.clear()
            params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
            result_request = requests.get(url=url_base, params=params)
            print(f"[Record no. {no}.]:{result_request.status_code}")

            blast_results.append(BlastResult(sam_record, "SUCCESS", rid, no, result_request.content.decode()))

            #with open(f'results_{rid}.txt', 'w') as save_file:
            #    save_file.write(result_request.content.decode())

            return True

        except Exception as e:
            print(f"Oops, something bad happened : {repr(e)}")
            return False











